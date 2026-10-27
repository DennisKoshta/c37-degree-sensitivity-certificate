#!/usr/bin/env python3
"""Translate a cert-degsens branch-and-bound certificate (certify.py) into VIPR 1.0 (exact MILP certificate format,
github.com/ambros-gleixner/VIPR), so the published checker viprchk can verify it independently. Stdlib only.

MILP: integer f_0..f_{N-1}; bounds f >= 0, f <= 1, f_0 <= 0, f_{e_i} >= 1; rows sum_{T<=S} (-1)^{|S|-|T|} f_T = 0
for |S| > d. A branch on f_j adds the assumptions f_j <= 0 and f_j >= 1 ({asm}); a leaf derives 0 >= delta > 0 as
a linear combination ({lin}) of rows and the current bounds; a branch node unsplits ({uns}) its two absurdities.

    python3 -I to_vipr.py IN.cert OUT.vipr
"""
import sys


def main(src, dst):
    lines = open(src).read().split("\n")
    if lines[-1] == "":
        lines.pop()
    _, _, n, d = lines[0].split()
    n, d = int(n), int(d)
    N = 1 << n
    pc = [bin(x).count("1") for x in range(N)]
    rows = [S for S in range(N) if pc[S] > d]
    out = open(dst, "w")
    w = out.write
    w("%% degree-sensitivity non-existence: n=%d, deg<=%d, f(0)=0, f(e_i)=1\nVER 1.0\n" % (n, d))
    w("VAR %d\n%s\n" % (N, " ".join(f"f{x}" for x in range(N))))
    w("INT %d\n%s\n" % (N, " ".join(str(x) for x in range(N))))
    w("OBJ min\n0\n")
    # bound constraints: lower bound of x at index x, upper bound at N + x
    nb = 2 * N
    w("CON %d %d\n" % (nb + len(rows), nb))
    lo_idx, hi_idx = {}, {}
    for x in range(N):
        lb = 1 if pc[x] == 1 else 0
        w(f"lb{x} G {lb} 1 {x} 1\n")
        lo_idx[x] = x
    for x in range(N):
        ub = 0 if x == 0 else 1
        w(f"ub{x} L {ub} 1 {x} 1\n")
        hi_idx[x] = N + x
    for r, S in enumerate(rows):
        terms = []
        T = S
        while True:
            terms.append(f"{T} {1 if (pc[S] - pc[T]) % 2 == 0 else -1}")
            if T == 0:
                break
            T = (T - 1) & S
        w(f"m{r} E 0 {len(terms)} {' '.join(terms)}\n")
    w("RTP infeas\nSOL 0\n")
    ders = []
    nxt = nb + len(rows)
    lo = [1 if pc[x] == 1 else 0 for x in range(N)]
    hi = [0 if x == 0 else 1 for x in range(N)]
    cur_lo, cur_hi = dict(lo_idx), dict(hi_idx)   # constraint index currently giving each bound

    def add(text):
        nonlocal nxt
        ders.append(text)
        nxt += 1
        return nxt - 1

    def leaf(tok):
        k = int(tok[1])
        c = {}
        lam = []
        for t in range(k):
            r, y = int(tok[2 + 2 * t]), int(tok[3 + 2 * t])
            lam.append((nb + r, y))
            S = rows[r]
            T = S
            while True:
                c[T] = c.get(T, 0) + (y if (pc[S] - pc[T]) % 2 == 0 else -y)
                if T == 0:
                    break
                T = (T - 1) & S
        # sum_r y_r (A_r f = 0) is c.f = 0; cancel c_j f_j with -c_j (f_j >= lo) if c_j < 0, -c_j (f_j <= hi) if
        # c_j > 0, leaving 0 >= -max_box c.f > 0
        rhs = 0
        for j, cj in sorted(c.items()):
            if cj < 0:
                lam.append((cur_lo[j], -cj))
                rhs += -cj * lo[j]
            elif cj > 0:
                lam.append((cur_hi[j], -cj))
                rhs += -cj * hi[j]
        assert rhs > 0, "leaf does not give an absurdity"
        body = " ".join(f"{i} {v}" for i, v in lam)
        return add(f"F{nxt} G {rhs} 0 {{ lin {len(lam)} {body} }} -1")

    pos = 1

    def node():
        nonlocal pos
        tok = lines[pos].split()
        pos += 1
        if tok[0] == "l":
            return leaf(tok)
        j = int(tok[1])
        saved = (lo[j], hi[j], cur_lo[j], cur_hi[j])
        a0 = add(f"A{nxt} L 0 1 {j} 1 {{ asm }} -1")
        hi[j], cur_hi[j] = 0, a0
        d0 = node()
        lo[j], hi[j], cur_lo[j], cur_hi[j] = saved
        a1 = add(f"A{nxt} G 1 1 {j} 1 {{ asm }} -1")
        lo[j], cur_lo[j] = 1, a1
        d1 = node()
        lo[j], hi[j], cur_lo[j], cur_hi[j] = saved
        return add(f"U{nxt} G 1 0 {{ uns {d0} {a0} {d1} {a1} }} -1")

    sys.setrecursionlimit(100000)
    node()
    assert pos == len(lines)
    w("DER %d\n" % len(ders))
    w("\n".join(ders) + "\n")
    out.close()
    print(f"wrote {dst}: {len(ders)} derivations")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
