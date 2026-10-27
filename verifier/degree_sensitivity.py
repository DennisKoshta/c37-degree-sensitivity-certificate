#!/usr/bin/env python3
"""Independent exact Möbius/sensitivity and interval verifier for C_37 (constant 37a).

Reference:
  Tao et al. optimizationproblems constant 37a (constants/37a.md):
    C_37 = sup_{n >= 2, deg(f) > 1} log(s(f)) / log(deg(f));
  Nisan, Szegedy (1994);
  Nisan, Wigderson (1995), Footnote 1 (Kushilevitz's 6-variable function);
  Buhrman, de Wolf (2002);
  Hatami, Kulkarni, Pankratov (2011), Example 5.4;
  Tal (2013); Huang (2019); Proskurin (2021).

File format (strictly 2 lines, UTF-8 without BOM, no CR or whitespace junk):
  Line 1: 'n d' (two canonical base-10 integers separated by a single space,
          with 2 <= n <= MAX_N and 2 <= d <= n; n is the number of variables /
          claimed sensitivity bound and d is the claimed degree bound).
  Line 2: lowercase hexadecimal string of length 2^n / 4 = 2^(n - 2).
          Bit x (for x = 0..2^n - 1) is bit (3 - (x % 4)) of hex digit (x // 4).
          f(x) = that bit in {0, 1}.
          Coordinate i (1-based, i = 1..n) of x is bit (i - 1) of x.

Mathematical definitions (confirmed against constants/37a.md):
  Multilinear polynomial representation over R:
      f(x_1, ..., x_n) = sum_{S subseteq [n]} c(S) prod_{i in S} x_i,
  where the exact integer Möbius coefficients are given by
      c(S) = sum_{T subseteq S} (-1)^{|S \\ T|} f(T).
  Real polynomial degree:
      deg(f) = max { |S| : c(S) != 0 }   (with deg(f) = 0 if c(S) = 0 for all S).
  Pointwise sensitivity at x in {0, 1}^n and maximum sensitivity:
      s(f)(x) = |{ i in {1, ..., n} : f(x) != f(x xor e_i) }|,
      s(f) = max_{x in {0, 1}^n} s(f)(x).
  Degree-sensitivity exponent (for deg(f) > 1):
      exponent(f) = log(s(f)) / log(deg(f)).
  Kushilevitz record comparison (default --beat="log_3 6"):
      s(f) > deg(f)^{log_3 6}  <=>  s(f)^{ln 3} > deg(f)^{ln 6}
                               <=>  ln(s(f)) * ln(3) > ln(deg(f)) * ln(6).
"""

import argparse
import decimal
from fractions import Fraction
import re
import sys
from typing import Any

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:  # pragma: no cover
    np = None  # type: ignore[assignment]
    HAS_NUMPY = False

try:
    from mpmath import iv as mp_iv

    HAS_MPMATH = True
except ImportError:  # pragma: no cover
    mp_iv = None  # type: ignore[assignment]
    HAS_MPMATH = False

MAX_N = 30
MAX_N_PURE_PYTHON = 14
PURE_PYTHON_CROSSCHECK_MAX_N = 14
INTERVAL_DPS = 100
OUTPUT_DIGITS = 25
MAX_STR_DIGITS = 2000
DEFAULT_BEAT = "log_3 6"

_CANONICAL_INT_RE = re.compile(r"\A-?(?:0|[1-9][0-9]*)\Z")
_LOWER_HEX_RE = re.compile(r"\A[0-9a-f]+\Z")
_BEAT_RE = re.compile(r"\A(?:0|[1-9][0-9]*)(?:\.[0-9]+|/[1-9][0-9]*)?\Z")
_LOG_BEAT_RE = re.compile(r"\Alog_([1-9][0-9]*)(?: ([1-9][0-9]*)|\(([1-9][0-9]*)\))\Z")

# Lookup table mapping each lowercase hex digit to its 4 {0, 1} truth-table entries
# for bits 3, 2, 1, 0: f(4*d + r) = (val >> (3 - r)) & 1.
_HEX_TO_01: dict[str, tuple[int, int, int, int]] = {
    f"{v:x}": (
        (v >> 3) & 1,
        (v >> 2) & 1,
        (v >> 1) & 1,
        v & 1,
    )
    for v in range(16)
}

if HAS_NUMPY:
    _POPCOUNT_16 = np.array([i.bit_count() for i in range(65536)], dtype=np.uint8)
else:  # pragma: no cover
    _POPCOUNT_16 = None  # type: ignore[assignment]


class _VerifierArgumentParser(argparse.ArgumentParser):
    def error(self, message: str):
        print(f"INVALID: {message}")
        raise SystemExit(1)


def parse_hex_to_01_pure_python(hex_str: str) -> list[int]:
    """Unpack lowercase hex string to a list of {0, 1} values in pure Python."""
    table = _HEX_TO_01
    return [val for ch in hex_str for val in table[ch]]


def parse_hex_to_01_numpy(hex_str: str, n: int) -> "np.ndarray":
    """Unpack lowercase hex string to a numpy uint8 array of {0, 1} values."""
    if not HAS_NUMPY:
        raise RuntimeError("numpy is not available")
    if n == 2:
        # 2^2 / 4 = 1 hex character (odd length, so bytes.fromhex is not applicable).
        return np.array(_HEX_TO_01[hex_str[0]], dtype=np.uint8)
    raw_bytes = bytes.fromhex(hex_str)
    return np.unpackbits(np.frombuffer(raw_bytes, dtype=np.uint8), bitorder="big")


def mobius_pure_python(f_list: list[int], n: int) -> list[int]:
    """Compute exact integer Möbius transform c(S) over Z in pure Python."""
    c = list(f_list)
    size = 1 << n
    if len(c) != size:
        raise ValueError(f"Expected length {size}, got {len(c)}")
    half = 1
    while half < size:
        step = half << 1
        for base in range(0, size, step):
            for j in range(base, base + half):
                c[j + half] -= c[j]
        half = step
    return c


def mobius_numpy(f_arr: "np.ndarray", n: int) -> "np.ndarray":
    """Compute exact integer Möbius transform c(S) over Z using numpy int64."""
    if not HAS_NUMPY:
        raise RuntimeError("numpy is not available")
    size = 1 << n
    if f_arr.shape != (size,):
        raise ValueError(f"Expected shape ({size},), got {f_arr.shape}")
    c = f_arr.astype(np.int64, copy=True)
    half = 1
    for _ in range(n):
        v = c.reshape(-1, 2, half)
        v[:, 1, :] -= v[:, 0, :]
        half <<= 1
    return c


def degree_from_mobius_pure_python(c_list: list[int]) -> int:
    """Return max |S| with c(S) != 0 (0 if all coefficients are 0)."""
    deg = 0
    for s, val in enumerate(c_list):
        if val != 0:
            bc = s.bit_count()
            if bc > deg:
                deg = bc
    return deg


def degree_from_mobius_numpy(c_arr: "np.ndarray", n: int) -> int:
    """Return max |S| with c(S) != 0 using 16-bit popcount blocks."""
    if not HAS_NUMPY:
        raise RuntimeError("numpy is not available")
    size = 1 << n
    if c_arr.shape != (size,):
        raise ValueError(f"Expected shape ({size},), got {c_arr.shape}")
    if n <= 16:
        nz = np.flatnonzero(c_arr)
        if nz.size == 0:
            return 0
        return int(_POPCOUNT_16[nz].max())

    c_2d = c_arr.reshape(1 << (n - 16), 65536)
    row_has_nz = np.any(c_2d != 0, axis=1)
    nz_rows = np.flatnonzero(row_has_nz)
    if nz_rows.size == 0:
        return 0
    best_deg = 0
    for hi in nz_rows:
        hi_int = int(hi)
        hi_bc = hi_int.bit_count()
        if hi_bc + 16 <= best_deg:
            continue
        nz_lo = np.flatnonzero(c_2d[hi_int])
        if nz_lo.size > 0:
            cand = hi_bc + int(_POPCOUNT_16[nz_lo].max())
            if cand > best_deg:
                best_deg = cand
    return best_deg


def sensitivity_pure_python(f_list: list[int], n: int) -> int:
    """Compute exact max sensitivity s(f) in pure Python."""
    size = 1 << n
    if len(f_list) != size:
        raise ValueError(f"Expected length {size}, got {len(f_list)}")
    best = 0
    for x in range(size):
        fx = f_list[x]
        sx = 0
        for i in range(n):
            if fx != f_list[x ^ (1 << i)]:
                sx += 1
        if sx > best:
            best = sx
    return best


def sensitivity_numpy(f_arr: "np.ndarray", n: int) -> int:
    """Compute exact max sensitivity s(f) using numpy uint8 coordinate slices."""
    if not HAS_NUMPY:
        raise RuntimeError("numpy is not available")
    size = 1 << n
    if f_arr.shape != (size,):
        raise ValueError(f"Expected shape ({size},), got {f_arr.shape}")
    f_u8 = f_arr.astype(np.uint8, copy=False)
    sens = np.zeros(size, dtype=np.uint8)
    half = 1
    for _ in range(n):
        fv = f_u8.reshape(-1, 2, half)
        sv = sens.reshape(-1, 2, half)
        diff = fv[:, 0, :] ^ fv[:, 1, :]
        sv[:, 0, :] += diff
        sv[:, 1, :] += diff
        half <<= 1
    return int(sens.max())


def _mpf_tuple_to_fraction(t: tuple[int, int, int, int]) -> Fraction:
    """Convert an internal mpmath (sign, man, exp, bc) tuple to an exact Fraction."""
    sign, man, exp, _bc = t
    if man == 0:
        return Fraction(0, 1)
    signed_man = -int(man) if sign else int(man)
    if exp >= 0:
        return Fraction(signed_man << int(exp), 1)
    return Fraction(signed_man, 1 << (-int(exp)))


def _exact_rational_log_ratio(s: int, deg: int) -> Fraction | None:
    """Return exact p/q if s^q == deg^p for integers p >= 0, q >= 1, else None."""
    if s < 1 or deg < 2:
        raise ValueError(f"Expected s >= 1 and deg >= 2, got s={s}, deg={deg}")
    if s == 1:
        return Fraction(0, 1)
    max_exp = max(s.bit_length(), deg.bit_length())
    for q in range(1, max_exp + 1):
        s_pow_q = s**q
        p = 1
        d_pow_p = deg
        while d_pow_p < s_pow_q:
            p += 1
            d_pow_p *= deg
        if d_pow_p == s_pow_q:
            return Fraction(p, q)
    return None


def _decimal_ln_exp_helpers(dps: int = INTERVAL_DPS):
    """Create directed-rounding Decimal contexts and rigorous ln/exp enclosures."""
    prec = max(int(dps), 80) + 15
    ctx_near = decimal.Context(prec=prec, rounding=decimal.ROUND_HALF_EVEN)
    ctx_lo = decimal.Context(prec=prec, rounding=decimal.ROUND_FLOOR)
    ctx_hi = decimal.Context(prec=prec, rounding=decimal.ROUND_CEILING)

    def ln_enclose(val_int: int) -> tuple[decimal.Decimal, decimal.Decimal]:
        if val_int == 1:
            z = decimal.Decimal(0)
            return z, z
        d = ctx_near.create_decimal(val_int)
        v = ctx_near.ln(d)
        ulp = v.next_plus(context=ctx_near) - v
        pad = ctx_hi.multiply(ulp, decimal.Decimal(4))
        return ctx_lo.subtract(v, pad), ctx_hi.add(v, pad)

    def exp_enclose(
        x_lo: decimal.Decimal, x_hi: decimal.Decimal
    ) -> tuple[decimal.Decimal, decimal.Decimal]:
        v_lo = ctx_near.exp(x_lo)
        ulp_lo = v_lo.next_plus(context=ctx_near) - v_lo
        pad_lo = ctx_hi.multiply(ulp_lo, decimal.Decimal(4))

        v_hi = ctx_near.exp(x_hi)
        ulp_hi = v_hi.next_plus(context=ctx_near) - v_hi
        pad_hi = ctx_hi.multiply(ulp_hi, decimal.Decimal(4))
        return ctx_lo.subtract(v_lo, pad_lo), ctx_hi.add(v_hi, pad_hi)

    return ctx_lo, ctx_hi, ln_enclose, exp_enclose


def _enclose_decimal_iv(
    s: int, deg: int, dps: int = INTERVAL_DPS
) -> tuple[Fraction, Fraction]:
    """Enclose log(s) / log(deg) using stdlib decimal directed-rounding intervals."""
    if s < 1 or deg < 2:
        raise ValueError(f"Expected s >= 1 and deg >= 2, got s={s}, deg={deg}")
    if s == 1:
        return Fraction(0, 1), Fraction(0, 1)

    ctx_lo, ctx_hi, ln_enclose, _ = _decimal_ln_exp_helpers(dps=dps)
    lns_lo, lns_hi = ln_enclose(s)
    lnd_lo, lnd_hi = ln_enclose(deg)

    r2_lo = Fraction(ctx_lo.divide(lns_lo, lnd_hi))
    r2_hi = Fraction(ctx_hi.divide(lns_hi, lnd_lo))

    e = 0
    m = s
    while m % deg == 0:
        m //= deg
        e += 1
    if m == 1:
        r1_lo = Fraction(e, 1)
        r1_hi = Fraction(e, 1)
    else:
        lnm_lo, lnm_hi = ln_enclose(m)
        rem_lo = ctx_lo.divide(lnm_lo, lnd_hi)
        rem_hi = ctx_hi.divide(lnm_hi, lnd_lo)
        r1_lo = Fraction(e, 1) + Fraction(rem_lo)
        r1_hi = Fraction(e, 1) + Fraction(rem_hi)

    lo = max(r1_lo, r2_lo)
    hi = min(r1_hi, r2_hi)
    exact = _exact_rational_log_ratio(s, deg)
    if exact is not None:
        if not (lo <= exact <= hi):
            raise AssertionError("Decimal enclosure failed to contain exact rational log ratio")
        lo, hi = exact, exact
    if lo > hi:
        raise AssertionError("Decimal exponent enclosures failed to overlap")
    return lo, hi


def _enclose_mpmath_iv(
    s: int, deg: int, dps: int = INTERVAL_DPS
) -> tuple[Fraction, Fraction]:
    """Enclose log(s) / log(deg) using mpmath.iv at >=80 digits in two formulations."""
    if s < 1 or deg < 2:
        raise ValueError(f"Expected s >= 1 and deg >= 2, got s={s}, deg={deg}")
    if s == 1:
        return Fraction(0, 1), Fraction(0, 1)

    old_dps = mp_iv.dps
    try:
        mp_iv.dps = max(int(dps), 80)
        ln_deg = mp_iv.log(mp_iv.mpf(deg))

        # Formulation 1: factor s = deg^e * m (deg not dividing m).
        e = 0
        m = s
        while m % deg == 0:
            m //= deg
            e += 1
        if m == 1:
            r_iv_1 = mp_iv.mpf(e)
        else:
            r_iv_1 = mp_iv.mpf(e) + mp_iv.log(mp_iv.mpf(m)) / ln_deg

        # Formulation 2: direct quotient log(s) / log(deg).
        r_iv_2 = mp_iv.log(mp_iv.mpf(s)) / ln_deg

        r1_lo = _mpf_tuple_to_fraction(r_iv_1._mpi_[0])
        r1_hi = _mpf_tuple_to_fraction(r_iv_1._mpi_[1])
        r2_lo = _mpf_tuple_to_fraction(r_iv_2._mpi_[0])
        r2_hi = _mpf_tuple_to_fraction(r_iv_2._mpi_[1])

        lo = max(r1_lo, r2_lo)
        hi = min(r1_hi, r2_hi)
        exact = _exact_rational_log_ratio(s, deg)
        if exact is not None:
            if not (lo <= exact <= hi):
                raise AssertionError(
                    "mpmath.iv enclosure failed to contain exact rational log ratio"
                )
            lo, hi = exact, exact
        if lo > hi:
            raise AssertionError("mpmath.iv exponent enclosures failed to overlap")
        return lo, hi
    finally:
        mp_iv.dps = old_dps


def enclose_exponent(
    s: int, deg: int, dps: int = INTERVAL_DPS
) -> tuple[Fraction, Fraction]:
    """Compute rigorous rational enclosure (exp_lo, exp_hi) of log(s) / log(deg)."""
    dec_lo, dec_hi = _enclose_decimal_iv(s, deg, dps=dps)
    if not HAS_MPMATH:
        return dec_lo, dec_hi

    mp_lo, mp_hi = _enclose_mpmath_iv(s, deg, dps=dps)
    if max(mp_lo, dec_lo) > min(mp_hi, dec_hi):
        raise AssertionError("mpmath.iv and Decimal exponent enclosures disagree")
    return mp_lo, mp_hi


def format_interval_outward(
    lo: Fraction, hi: Fraction, digits: int = OUTPUT_DIGITS
) -> tuple[str, str, Fraction, Fraction]:
    """Format rational interval [lo, hi] (>= 0) with outward rounding to `digits` decimal places."""
    if lo < 0 or hi < lo:
        raise ValueError(f"Expected 0 <= lo <= hi, got [{lo}, {hi}]")
    scale = 10**digits
    scaled_lo = (lo.numerator * scale) // lo.denominator
    scaled_hi = (hi.numerator * scale + hi.denominator - 1) // hi.denominator

    lo_int, lo_frac_part = divmod(scaled_lo, scale)
    hi_int, hi_frac_part = divmod(scaled_hi, scale)
    lo_str = f"{lo_int}.{lo_frac_part:0{digits}d}"
    hi_str = f"{hi_int}.{hi_frac_part:0{digits}d}"
    lo_dec_frac = Fraction(scaled_lo, scale)
    hi_dec_frac = Fraction(scaled_hi, scale)
    return lo_str, hi_str, lo_dec_frac, hi_dec_frac


def _exact_same_log_ratio(s: int, deg: int, arg_a: int, base_b: int) -> bool:
    """Return True if log(s)/log(deg) == log(arg_a)/log(base_b) by exact integer powers."""
    r1 = _exact_rational_log_ratio(s, deg)
    r2 = _exact_rational_log_ratio(arg_a, base_b)
    if r1 is not None or r2 is not None:
        return r1 == r2
    max_exp = max(s.bit_length(), deg.bit_length(), arg_a.bit_length(), base_b.bit_length())
    for k2 in range(1, max_exp + 1):
        s_k2 = s**k2
        d_k2 = deg**k2
        a_k1 = arg_a
        b_k1 = base_b
        while a_k1 < s_k2 and b_k1 < d_k2:
            a_k1 *= arg_a
            b_k1 *= base_b
        if a_k1 == s_k2 and b_k1 == d_k2:
            return True
    return False


def _certify_beats_log_threshold(
    s: int, deg: int, base_b: int, arg_a: int, dps: int = INTERVAL_DPS
) -> bool:
    """Certify s > deg^{log_{base_b}(arg_a)} via rigorous power and log enclosures."""
    if s <= 1:
        return False

    exact_eq = _exact_same_log_ratio(s, deg, arg_a, base_b)
    cur_dps = max(int(dps), 80)

    while True:
        ctx_lo, ctx_hi, ln_enclose, exp_enclose = _decimal_ln_exp_helpers(dps=cur_dps)
        lns_lo, lns_hi = ln_enclose(s)
        lnd_lo, lnd_hi = ln_enclose(deg)
        lna_lo, lna_hi = ln_enclose(arg_a)
        lnb_lo, lnb_hi = ln_enclose(base_b)

        # Log-product comparison: ln(s) * ln(b) > ln(deg) * ln(a).
        lhs_log_lo = ctx_lo.multiply(lns_lo, lnb_lo)
        lhs_log_hi = ctx_hi.multiply(lns_hi, lnb_hi)
        rhs_log_lo = ctx_lo.multiply(lnd_lo, lna_lo)
        rhs_log_hi = ctx_hi.multiply(lnd_hi, lna_hi)
        dec_log_beats = lhs_log_lo > rhs_log_hi
        dec_log_below = lhs_log_hi < rhs_log_lo

        # Power comparison: s^{ln b} > deg^{ln a}.
        s_pow_lo, s_pow_hi = exp_enclose(lhs_log_lo, lhs_log_hi)
        d_pow_lo, d_pow_hi = exp_enclose(rhs_log_lo, rhs_log_hi)
        dec_pow_beats = s_pow_lo > d_pow_hi
        dec_pow_below = s_pow_hi < d_pow_lo

        exp_lo, exp_hi = _enclose_decimal_iv(s, deg, dps=cur_dps)
        beat_lo, beat_hi = _enclose_decimal_iv(arg_a, base_b, dps=cur_dps)
        dec_ratio_beats = exp_lo > beat_hi
        dec_ratio_below = exp_hi < beat_lo

        dec_beats = dec_log_beats and dec_pow_beats and dec_ratio_beats
        dec_separated = dec_beats or (
            dec_log_below and dec_pow_below and dec_ratio_below
        )

        if not HAS_MPMATH:
            if exact_eq:
                if dec_beats:
                    raise AssertionError("Decimal beat check certified True on exact equality")
                return False
            if dec_separated or cur_dps >= 3200:
                return dec_beats
            cur_dps *= 2
            continue

        old_dps = mp_iv.dps
        try:
            mp_iv.dps = cur_dps
            mp_s = mp_iv.mpf(s)
            mp_d = mp_iv.mpf(deg)
            mp_a = mp_iv.mpf(arg_a)
            mp_b = mp_iv.mpf(base_b)

            ln_s = mp_iv.log(mp_s)
            ln_d = mp_iv.log(mp_d)
            ln_a = mp_iv.log(mp_a)
            ln_b = mp_iv.log(mp_b)

            # Direct power test: s^(ln b) > deg^(ln a) in mpmath.iv.
            s_pow_iv = mp_s ** ln_b
            d_pow_iv = mp_d ** ln_a
            mp_sp_lo = _mpf_tuple_to_fraction(s_pow_iv._mpi_[0])
            mp_sp_hi = _mpf_tuple_to_fraction(s_pow_iv._mpi_[1])
            mp_dp_lo = _mpf_tuple_to_fraction(d_pow_iv._mpi_[0])
            mp_dp_hi = _mpf_tuple_to_fraction(d_pow_iv._mpi_[1])
            mp_pow_beats = mp_sp_lo > mp_dp_hi
            mp_pow_below = mp_sp_hi < mp_dp_lo

            # Cross-multiplied log enclosure test: ln(s)*ln(b) - ln(deg)*ln(a) > 0.
            log_diff_iv = ln_s * ln_b - ln_d * ln_a
            mp_ld_lo = _mpf_tuple_to_fraction(log_diff_iv._mpi_[0])
            mp_ld_hi = _mpf_tuple_to_fraction(log_diff_iv._mpi_[1])
            mp_log_beats = mp_ld_lo > 0
            mp_log_below = mp_ld_hi < 0

            # Exponent ratio enclosure test: (ln s / ln deg)_lo > (ln a / ln b)_hi.
            mp_exp_lo, mp_exp_hi = _enclose_mpmath_iv(s, deg, dps=cur_dps)
            mp_beat_lo, mp_beat_hi = _enclose_mpmath_iv(arg_a, base_b, dps=cur_dps)
            mp_ratio_beats = mp_exp_lo > mp_beat_hi
            mp_ratio_below = mp_exp_hi < mp_beat_lo

            mp_beats = mp_pow_beats and mp_log_beats and mp_ratio_beats
            mp_separated = mp_beats or (
                mp_pow_below and mp_log_below and mp_ratio_below
            )
        finally:
            mp_iv.dps = old_dps

        if exact_eq:
            if mp_beats or dec_beats:
                raise AssertionError("Beat check certified True on exact equality")
            return False
        if (mp_separated and dec_separated) or cur_dps >= 3200:
            if mp_beats != dec_beats:
                raise AssertionError("mpmath.iv and Decimal beat certifications disagree")
            return mp_beats
        cur_dps *= 2


def _certify_beats_fraction_threshold(
    s: int,
    deg: int,
    beat_frac: Fraction,
    exp_lo: Fraction | None = None,
    exp_hi: Fraction | None = None,
    dps: int = INTERVAL_DPS,
) -> tuple[bool, Fraction, Fraction]:
    """Certify log(s)/log(deg) > beat_frac from rigorous enclosure, increasing dps if needed."""
    cur_dps = max(int(dps), 80)
    if exp_lo is None or exp_hi is None:
        exp_lo, exp_hi = enclose_exponent(s, deg, dps=cur_dps)

    while True:
        if exp_lo > beat_frac:
            return True, exp_lo, exp_hi
        if exp_hi <= beat_frac:
            return False, exp_lo, exp_hi
        if cur_dps >= 3200:
            return False, exp_lo, exp_hi
        cur_dps *= 2
        exp_lo, exp_hi = enclose_exponent(s, deg, dps=cur_dps)


def check_beats_threshold(
    s: int,
    deg: int,
    exp_lo: Fraction | None = None,
    lo_dec_frac: Fraction | None = None,
    beat: str | Fraction | None = DEFAULT_BEAT,
    dps: int = INTERVAL_DPS,
    exp_hi: Fraction | None = None,
) -> tuple[bool, str]:
    """Validate `--beat` and return (beats, canonical_beat_str)."""
    if s < 1 or deg < 2:
        raise ValueError(f"Expected s >= 1 and deg >= 2, got s={s}, deg={deg}")

    if beat is None:
        beat = DEFAULT_BEAT

    if isinstance(beat, Fraction):
        if beat < 0:
            raise ValueError("beat threshold must be non-negative")
        beats, _, _ = _certify_beats_fraction_threshold(
            s, deg, beat, exp_lo=exp_lo, exp_hi=exp_hi, dps=dps
        )
        return beats, str(beat)

    if not isinstance(beat, str):
        raise ValueError("beat threshold must be a string or Fraction")
    if not beat or beat != beat.strip():
        raise ValueError("whitespace junk in --beat value")
    if len(beat) > MAX_STR_DIGITS:
        raise ValueError("--beat value exceeds digit limit")

    log_match = _LOG_BEAT_RE.match(beat)
    if log_match:
        b_str = log_match.group(1)
        a_str = log_match.group(2) or log_match.group(3)
        if len(b_str) > 10 or len(a_str) > 10:
            raise ValueError("log base or argument in --beat is too large")
        base_b = int(b_str)
        arg_a = int(a_str)
        if base_b < 2 or arg_a < 2:
            raise ValueError("require base >= 2 and argument >= 2 in log_b a")
        beats = _certify_beats_log_threshold(s, deg, base_b, arg_a, dps=dps)
        return beats, beat

    if any(c.isspace() for c in beat):
        raise ValueError("whitespace junk in --beat value")
    if not _BEAT_RE.match(beat):
        raise ValueError(
            "expected 'log_b a' or canonical unsigned decimal / rational 'p/q' for --beat"
        )
    beat_frac = Fraction(beat)
    init_dps = max(int(dps), len(beat) + 50, 80)
    beats, _, _ = _certify_beats_fraction_threshold(
        s,
        deg,
        beat_frac,
        exp_lo=exp_lo if init_dps <= dps else None,
        exp_hi=exp_hi if init_dps <= dps else None,
        dps=init_dps,
    )
    return beats, beat


def verify_degree_sensitivity_lines(
    lines: list[str],
    beat: str | Fraction | None = DEFAULT_BEAT,
    max_n: int = MAX_N,
    use_numpy: bool | None = None,
) -> tuple[bool, dict[str, Any] | None, str]:
    """Parse and verify a 2-line C_37 Boolean function construction, returning (ok, details, msg)."""
    if len(lines) != 2:
        return (
            False,
            None,
            f"INVALID: expected exactly 2 lines ('n d' and hex string), got {len(lines)} lines",
        )

    header_line, hex_line = lines[0], lines[1]

    if "\ufeff" in header_line or "\ufeff" in hex_line:
        return False, None, "INVALID: UTF-8 BOM junk in input"
    if "\r" in header_line or "\r" in hex_line:
        return False, None, "INVALID: carriage return whitespace junk in input"

    if not header_line:
        return False, None, "INVALID: line 1: empty line for 'n d'"
    if (
        header_line != header_line.strip()
        or any(c.isspace() and c != " " for c in header_line)
        or "  " in header_line
    ):
        return False, None, "INVALID: line 1: whitespace junk in 'n d'"

    header_tokens = header_line.split(" ")
    if len(header_tokens) != 2:
        return (
            False,
            None,
            f"INVALID: line 1: expected 'n d' (two canonical integers), got {len(header_tokens)} tokens",
        )

    n_str, d_str = header_tokens[0], header_tokens[1]
    if not _CANONICAL_INT_RE.match(n_str):
        return (
            False,
            None,
            f"INVALID: line 1: expected canonical base-10 integer n, got {n_str!r}",
        )
    if not _CANONICAL_INT_RE.match(d_str):
        return (
            False,
            None,
            f"INVALID: line 1: expected canonical base-10 integer d, got {d_str!r}",
        )
    if len(n_str) > 10:
        return (
            False,
            None,
            f"INVALID: line 1: n exceeds limit ({max_n}); too large for this verifier",
        )
    if len(d_str) > 10:
        return (
            False,
            None,
            f"INVALID: line 1: d exceeds limit ({max_n}); too large for this verifier",
        )

    n = int(n_str)
    d = int(d_str)
    if n < 2:
        return (
            False,
            None,
            f"INVALID: line 1: require n >= 2 so 2^n/4 is an integer, got n={n}",
        )
    if n > max_n:
        return (
            False,
            None,
            f"INVALID: line 1: n = {n} exceeds limit ({max_n}); too large for this verifier",
        )
    if d < 0 or d > n:
        return (
            False,
            None,
            f"INVALID: line 1: claimed degree d={d} out of range 2..{n}",
        )

    if use_numpy is None:
        use_numpy_flag = HAS_NUMPY
    else:
        use_numpy_flag = bool(use_numpy)
        if use_numpy_flag and not HAS_NUMPY:
            return False, None, "INVALID: numpy requested but not available"

    if not use_numpy_flag and n > MAX_N_PURE_PYTHON:
        return (
            False,
            None,
            f"INVALID: line 1: n = {n} exceeds pure-Python limit ({MAX_N_PURE_PYTHON})",
        )

    if not hex_line:
        return False, None, "INVALID: line 2: empty hex string"
    if any(c.isspace() for c in hex_line):
        return False, None, "INVALID: line 2: whitespace junk in hex string"
    if not _LOWER_HEX_RE.match(hex_line):
        return (
            False,
            None,
            "INVALID: line 2: hex string must contain only lowercase hex digits 0-9a-f",
        )

    expected_hex_len = 1 << (n - 2)
    if len(hex_line) != expected_hex_len:
        return (
            False,
            None,
            f"INVALID: line 2: expected hex string of length {expected_hex_len} (2^{n}/4), "
            f"got {len(hex_line)}",
        )

    try:
        if use_numpy_flag:
            f_arr = parse_hex_to_01_numpy(hex_line, n)
            c_arr = mobius_numpy(f_arr, n)

            # Verify Möbius identities c(empty) = f(0) and sum_S c(S) = f(1, ..., 1).
            if int(c_arr[0]) != int(f_arr[0]):
                raise AssertionError(
                    f"Möbius constant term mismatch: c(0) = {int(c_arr[0])} != f(0) = {int(f_arr[0])}"
                )
            if int(c_arr.sum()) != int(f_arr[-1]):
                raise AssertionError(
                    f"Möbius sum identity failed: sum_S c(S) = {int(c_arr.sum())} != f(2^n-1) = {int(f_arr[-1])}"
                )

            deg = degree_from_mobius_numpy(c_arr, n)
            s = sensitivity_numpy(f_arr, n)

            # Cross-check numpy vs pure-Python for n <= 14.
            if n <= PURE_PYTHON_CROSSCHECK_MAX_N:
                f_py = parse_hex_to_01_pure_python(hex_line)
                if f_arr.tolist() != f_py:
                    raise AssertionError("numpy and pure-Python hex unpacking disagreed")
                c_py = mobius_pure_python(f_py, n)
                if c_arr.tolist() != c_py:
                    raise AssertionError("numpy and pure-Python Möbius transform disagreed")
                deg_py = degree_from_mobius_pure_python(c_py)
                if deg != deg_py:
                    raise AssertionError("numpy and pure-Python degree disagreed")
                s_py = sensitivity_pure_python(f_py, n)
                if s != s_py:
                    raise AssertionError("numpy and pure-Python sensitivity disagreed")
        else:
            f_py = parse_hex_to_01_pure_python(hex_line)
            c_py = mobius_pure_python(f_py, n)
            if c_py[0] != f_py[0]:
                raise AssertionError(
                    f"Möbius constant term mismatch: c(0) = {c_py[0]} != f(0) = {f_py[0]}"
                )
            if sum(c_py) != f_py[-1]:
                raise AssertionError(
                    f"Möbius sum identity failed: sum_S c(S) = {sum(c_py)} != f(2^n-1) = {f_py[-1]}"
                )
            deg = degree_from_mobius_pure_python(c_py)
            s = sensitivity_pure_python(f_py, n)
    except MemoryError as exc:  # pragma: no cover
        return False, None, f"INVALID: out of memory verifying n={n}: {exc}"

    if deg <= 1:
        return (
            False,
            None,
            f"INVALID: degenerate function with deg={deg} <= 1 (require deg >= 2)",
        )
    if d <= 1:
        return (
            False,
            None,
            f"INVALID: line 1: claimed degree d={d} <= 1 (require d >= 2)",
        )
    if deg < d:
        return (
            False,
            None,
            f"INVALID: function has lower degree deg={deg} than claimed d={d}",
        )
    if deg > d:
        return (
            False,
            None,
            f"INVALID: function has higher degree deg={deg} than claimed d={d}",
        )
    if s < 1:
        return (
            False,
            None,
            f"INVALID: degenerate function with s={s} < 1",
        )

    exp_lo, exp_hi = enclose_exponent(s, deg)
    lo_str, hi_str, lo_dec_frac, hi_dec_frac = format_interval_outward(
        exp_lo, exp_hi, digits=OUTPUT_DIGITS
    )

    try:
        beats, beat_str = check_beats_threshold(
            s, deg, exp_lo=exp_lo, lo_dec_frac=lo_dec_frac, beat=beat, exp_hi=exp_hi
        )
    except (ValueError, ZeroDivisionError) as exc:
        return False, None, f"INVALID: invalid --beat threshold {beat!r}: {exc}"

    if beats and beat is not None and not (
        isinstance(beat, str) and _LOG_BEAT_RE.match(beat)
    ):
        beat_frac = beat if isinstance(beat, Fraction) else Fraction(beat)
        if exp_lo <= beat_frac:
            init_dps = max(INTERVAL_DPS, len(str(beat)) + 50, 80)
            _, exp_lo, exp_hi = _certify_beats_fraction_threshold(
                s, deg, beat_frac, dps=init_dps
            )
            lo_str, hi_str, lo_dec_frac, hi_dec_frac = format_interval_outward(
                exp_lo, exp_hi, digits=OUTPUT_DIGITS
            )
        if lo_dec_frac <= beat_frac:
            diff = exp_lo - beat_frac
            digits = max(
                OUTPUT_DIGITS + 1,
                len(str(diff.denominator // diff.numerator)),
            )
            while True:
                lo_str, hi_str, lo_dec_frac, hi_dec_frac = format_interval_outward(
                    exp_lo, exp_hi, digits=digits
                )
                if lo_dec_frac > beat_frac:
                    break
                digits += 1

    msg = f"VALID n_vars={n} deg={deg} s={s} exponent_lo={lo_str}"
    if beats:
        msg += f" BEATS {beat_str}"

    details: dict[str, Any] = {
        "n": n,
        "n_vars": n,
        "claimed_d": d,
        "deg": deg,
        "s": s,
        "exp_lo": exp_lo,
        "exp_hi": exp_hi,
        "lo_str": lo_str,
        "hi_str": hi_str,
        "lo_dec": lo_dec_frac,
        "hi_dec": hi_dec_frac,
        "beats": beats,
        "beat_str": beat_str,
    }
    return True, details, msg


def verify_degree_sensitivity_file(
    filepath: str,
    beat: str | Fraction | None = DEFAULT_BEAT,
    max_n: int = MAX_N,
    use_numpy: bool | None = None,
) -> tuple[bool, dict[str, Any] | None, str]:
    """Read a C_37 degree-sensitivity construction file and verify it."""
    try:
        with open(filepath, "r", encoding="utf-8", newline="") as f:
            raw = f.read()
    except (OSError, UnicodeDecodeError) as exc:
        return False, None, f"INVALID: cannot read file: {exc}"

    if not raw:
        return False, None, "INVALID: empty file"
    if "\ufeff" in raw:
        return False, None, "INVALID: UTF-8 BOM junk in file"
    if "\r" in raw:
        return False, None, "INVALID: carriage return whitespace junk in file"

    # Allow at most one trailing newline at EOF; reject trailing blank lines.
    body = raw[:-1] if raw.endswith("\n") else raw
    if not body:
        return False, None, "INVALID: empty file"
    lines = body.split("\n")
    return verify_degree_sensitivity_lines(
        lines, beat=beat, max_n=max_n, use_numpy=use_numpy
    )


def main(argv: list[str] | None = None) -> int:
    parser = _VerifierArgumentParser(
        description="Verify a C_37 (constant 37a) degree-sensitivity Boolean function file."
    )
    parser.add_argument("file", help="Path to the 2-line C_37 construction file")
    parser.add_argument(
        "--beat",
        type=str,
        default=DEFAULT_BEAT,
        help="Threshold to beat (default: 'log_3 6'); prints BEATS only if strictly exceeded",
    )
    args = parser.parse_args(argv)

    ok, _details, message = verify_degree_sensitivity_file(args.file, beat=args.beat)
    print(message)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
