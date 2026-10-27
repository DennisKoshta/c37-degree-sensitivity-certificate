# No Boolean function of degree 4 has sensitivity 10

*A computer-assisted non-existence result related to the degree–sensitivity exponent (constant 37a in Tao et al.'s
optimization-constants table). Version of 2026-10-26 (first draft 2026-10-23). Not peer reviewed; see the AI-use
disclosure at the end.*

## 1. Statement

For $f:\{0,1\}^n\to\{0,1\}$ let $\deg f$ be the degree of its unique multilinear real polynomial. Let
$s(f)=\max_x |\{i : f(x)\ne f(x\oplus e_i)\}|$ be its (maximum) sensitivity.

**Theorem.** Every Boolean function $f$ with $\deg f\le 4$ has $s(f)\le 9$. Equivalently, there is no
$g:\{0,1\}^{10}\to\{0,1\}$ with $\deg g\le 4$, $g(0)=0$ and $g(e_i)=1$ for $i=1,\dots,10$.

The bound is attained. $\mathrm{NAE}_3\circ\mathrm{NAE}_3$ on 9 variables, where
$\mathrm{NAE}_3(x)=x_1+x_2+x_3-x_1x_2-x_1x_3-x_2x_3$, has degree $4$ and sensitivity $9$. Together with the
same method for smaller degrees (§4.4), the largest sensitivity of a Boolean function of degree at most $d$ is

$$ m(2)=3,\qquad m(3)=6,\qquad m(4)=9 . $$

The value $m(3)=6$ is Kushilevitz's function, and Maitra et al. [MMST] ruled out 7 earlier, with an uncertified ILP.

## 2. Context

**The constant.** Tao et al.'s table entry 37a [T37] defines $C_{37}=\sup \log s(f)/\log\deg f$ over Boolean $f$
with $\deg f>1$. It records $\log_3 6\approx 1.63093\le C_{37}\le 2$.
- The lower bound is Kushilevitz's 6-variable function: degree 3, sensitivity 6. It appears in footnote 1 of
  Nisan–Wigderson (1995) [NW95] and in Hatami–Kulkarni–Pankratov (2011, Example 5.4) [HKP11].
- The upper bound $2$ follows from $s\le bs$ and $bs(f)=O(\deg(f)^2)$, by Nisan–Szegedy (1994) [NS94] and Tal (2013)
  [Tal13]. Proskurin (2021) [P21] and Wellens (2022) [W22] improved the constant factor.

On 2026-10-23 the entry still lists $\log_3 6$ (repository `teorth/optimizationproblems` at commit `46dbf03`,
2026-10-17). A search of that repository's pull requests and issues found none that changes 37a.

**Why $(n,d)=(10,4)$.** Restrict to the sensitive coordinates of a maximally sensitive point and XOR-shift (Lemma 1
below). Any $f$ with $s(f)=n$ and $\deg f=d$ then yields a fully sensitive $g$ on $n$ variables with $g(0)=0$,
$g(e_i)=1$ and $\deg g\le d$. Composing $g$ with itself $k$ times keeps full sensitivity at $0$, with $s=n^k$ and
$\deg\le d^k$, so $C_{37}\ge \log n/\log d$ [HKP11, Ex. 5.4; MMST, Thm 13].

So a single fully sensitive function with $n>d^{\log_3 6}$ would beat the record. For $d=4$ that means $n\ge 10$,
which would give $\log 10/\log 4\approx 1.661$. A symmetrisation (linear-programming) argument rules out $n\ge 11$ for
$d=4$ [W22, Table 1], so $(10,4)$ was the only open case at degree 4.

**Prior work on $(10,4)$.**
- Proskurin [P21, §5; arXiv 2021-01] singled out this function. He proved (his Theorem 11) that if it exists, its
  symmetrisation is unique:
  $p(k)=0,1,1,\tfrac{7}{12},\tfrac16,0,\tfrac16,\tfrac{7}{12},1,1,0$ for $|x|=k=0,\dots,10$. He left existence open.
- Maitra, Mukherjee, Stănică and Tang [MMST; arXiv 2107.11205, 2021-07; Discrete Appl. Math. 365 (2025) 138–159]
  showed there is no *rotation-symmetric* such function (their "$[10,1,5]$"). They also showed by Gurobi ILP, without a
  certificate, that no fully sensitive function exists on 7 variables with degree 3.

As far as we can tell, the unrestricted $(10,4)$ case had not been settled.

**Re-check (2026-10-26).** `constants/37a.md` is unchanged at `teorth/optimizationproblems` main `e95af44`
(last edit to the file 2026-03-09). No open pull request touches it, and searching the repository's issues and pull
requests for "37a" and "sensitivity" finds nothing. An arXiv search for submissions from 2026-10-20 to 2026-10-27
(sensitivity and degree, Kushilevitz, Boolean function degree) found nothing related.

**Literature check (2026-10-23).**
- We searched arXiv (2025–2026) and the web for sensitivity versus degree, Kushilevitz, Proskurin's function and
  "degree 4 / sensitivity 10", and looked at the 37a page and the Epoch AI open-problem page.
- Nothing found settles $(10,4)$ or improves $\log_3 6$.
- Balodis (arXiv 2103.05593, 2021) improves the Kushilevitz–Nisan *certificate complexity* versus degree separation to
  $C(f)\ge\deg(f)^{2-o(1)}$, and $C$ versus $s$ to exponent $3-o(1)$. These do not bear on $s$ versus $\deg$.
- Li and Zhang (arXiv 2610.17006 and 2610.17100, 2026-10) compute sensitivities of rotation-symmetric functions but do
  not address this question.
- A search snippet mentions a talk by Smotrovs and Zaicevs at Theory Days 2026 (April) on low-degree representing
  polynomials of high-sensitivity functions. We could not open the slides (HTTP 403), so we cannot rule out overlap
  with this note.

## 3. Proof

**Lemma 1 (reduction).** If $\deg f\le d$ and $s(f)\ge n$, there is $g:\{0,1\}^n\to\{0,1\}$ with $\deg g\le d$,
$g(0)=0$ and $g(e_i)=1$ for all $i$.

*Proof.* Take $x$ with $s(f)(x)\ge n$ and $n$ coordinates $S$ sensitive at $x$. Fix the other coordinates to their
values in $x$, which does not raise the degree. Replace $y_i$ by $1-y_i$ where $x_i=1$, so that $x$ maps to $0$; this
preserves degree. If the resulting value at $0$ is $1$, replace the function by $1-g$. $\square$

Restricting such a $g$ on $n+1$ variables to $x_{n+1}=0$ gives one on $n$ variables. So the theorem is equivalent to
the non-existence of $g$ with $n=10$, $d=4$, and this also excludes every $n\ge 10$.

**Integer program.** Write $g_T$ for $g(1_T)$, $T\subseteq[10]$. The Möbius coefficients are
$a_S=\sum_{T\subseteq S}(-1)^{|S|-|T|}g_T$, and $\deg g\le 4$ iff $a_S=0$ for all $|S|\ge 5$. So $g$ exists iff the
system

$$ A g = 0,\qquad g\in\{0,1\}^{1024},\qquad g_\emptyset=0,\qquad g_{\{i\}}=1\ (1\le i\le 10) \tag{$*$}$$

is solvable. Here $A$ has one row per $S$ with $|S|\ge5$ (638 rows), with entries $(-1)^{|S|-|T|}$ for $T\subseteq S$.

The LP relaxation of $(*)$ is *feasible*: $g_T=p(|T|)$ with Proskurin's $p$ is a solution. So a pure LP-duality
(symmetrisation) proof cannot work, and integrality has to be used. Two cheaper obstructions also fail:
- reducing mod 2 gives a system consistent with the Reed–Muller code $RM(4,10)$ (rank 576 on 660 unknowns);
- the Hamming weight forced by $p$ is $320\equiv 0 \pmod{64}$, so a weight-divisibility argument doesn't apply.

**Certificate.** We give a binary branch-and-bound tree over the variables $g_T$. Each leaf $\ell$ comes with an
integer vector $y^{(\ell)}\in\mathbb{Z}^{638}$. Let $B_\ell$ be the box $\{0\le g\le 1\}$ intersected with the root
fixings and the branch fixings on the path to $\ell$. The leaf condition is

$$ \max_{g\in B_\ell}\ \big(y^{(\ell)\top}A\big)\,g \;<\;0 . $$

The left side is a linear function maximised over a box, so it is computed exactly coordinate by coordinate. Any $g$
with $Ag=0$ has $(y^\top A)g=0$, so no solution of $(*)$ lies in $B_\ell$. The leaf boxes cover every 0/1 point
satisfying the root fixings, so $(*)$ has no solution. $\square$

The tree comes from plain LP branch-and-bound: HiGHS 1.15.1 LP solves, most-fractional branching, no cutting planes,
no symmetry reduction. At each infeasible node, the dual of a phase-1 LP is scaled by a power of two (at most
$2^{32}$), rounded to integers and re-checked exactly. The certificate has:
- 957 nodes, 478 branchings, 479 leaves, depth 15;
- 638 nonzero multipliers at every leaf, with $\max|y|=2^{24}$.

## 4. How it was checked

1. **Our checker.** `checker/check_cert.py` uses only the Python standard library and exact integers. It rebuilds $A$ from
   $(n,d)=(10,4)$ alone, walks the tree in preorder, and checks every leaf inequality. It reports
   `VALID ... (479 leaves, depth 15, sha256 7dc38b37…)` in 15 s. It rejects a copy with one leaf negated and a copy
   with one leaf deleted. This checker was written by the same agent that produced the certificate.
2. **Independent proof checker.** `checker/to_vipr.py` translates the tree into the VIPR 1.0 format for exact MILP
   certificates [VIPR]:
   - bounds and Möbius rows as constraints;
   - each branch as two `asm` assumptions closed by `uns`;
   - each leaf as a `lin` combination of rows and bounds giving $0\ge\delta>0$.

   The published checker `viprchk` reports **"Infeasibility verified."** in 12.0 s of CPU. We built it from
   `github.com/ambros-gleixner/VIPR` at commit `f221c3fd2274e8e859e8e24505278d6f7c1bb8ff` against GMP 6.3.0. It
   rejects a copy with one leaf multiplier negated ("Failed to derive constraint").

   `viprchk` checks the VIPR file against the problem *as stated in that file*. That problem section (1024 integer
   variables, 2048 bounds, 638 equality rows) was produced by our translator, so a reader should satisfy themselves
   that it is $(*)$. It is short and regular enough to audit directly.
3. **Re-runs.**
   - On 2026-10-23, from the committed compressed certificate: decompression gives the same sha256, and both checks
     pass again (checker VALID in 15 s; `viprchk` "Infeasibility verified." in 12.0 s). The regenerated VIPR file has
     the same sha256 as on 2026-10-22.
   - On 2026-10-26, the commands in `README.md` were run from a scratch copy of this repository, with `viprchk`
     rebuilt from source. The results were the same. Regenerating the certificate with `generator/certify.py` gave a
     byte-identical file (sha256 `7dc38b37…`).
4. **Sharpness and the smaller degrees.**
   - The independent C37 verifier `verifier/degree_sensitivity.py` (written by a separate AI agent from the table's
     definition, without the search code) confirms that $\mathrm{NAE}_3\circ\mathrm{NAE}_3$ has degree $4$ and
     sensitivity $9$ (`VALID n_vars=9 deg=4 s=9`). That verifier checks explicit functions; it does not check
     non-existence certificates.
   - The same pipeline certifies $(n,d)=(4,2)$ and $(7,3)$ infeasible, each with a single LP-infeasible leaf, and
     both pass `checker/check_cert.py` and `viprchk`. With NAE ($s=3$, degree 2) and Kushilevitz ($s=6$, degree 3), this gives
     $m(2)=3$ and $m(3)=6$.
   - $(10,4)$ is the first case in this family where the LP relaxation is feasible.

## 5. What is not claimed

- **No record.** Constant 37a stays at $\log_3 6$. The theorem removes the one rung, $(10,4)$, that could have beaten
  it at degree 4. It says nothing about larger degrees: a fully sensitive function with $(n,d)=(14,5)$ or $(15,5)$,
  $(19..21,6)$, $(24..29,7)$, … would still beat $\log_3 6$, and none of these is excluded.
- **Search results are not proofs.** We also searched for fully sensitive functions invariant under permutation groups
  at those larger rungs. Those logs are not included in this repository. Every completed search was infeasible according to HiGHS, but
  those verdicts carry no certificates, many cases timed out, and none of it is part of the theorem.
- **Novelty is not established.** We found no prior proof of the unrestricted $(10,4)$ case (§2), but we could not
  read one possibly related 2026 talk.
- **No human has checked the proof.** The only checks are the computer checks above.

## 6. Reproduction

See `README.md` for the commands, which have been run from a scratch copy, and `SHA256SUMS` for the file hashes.

The certificate format is documented at the top of `generator/certify.py`. In it, `b M` branches on $g_M$ ($g_M=0$
subtree first), and `l k r_1 y_1 … r_k y_k` is a leaf, where row $r$ is the $r$-th set $S$ with $|S|\ge5$ in
increasing bitmask order. Bit $i-1$ of a mask is coordinate $i$.

| File | sha256 |
| --- | --- |
| `certificates/d10_4.cert.xz` (722,080 bytes) | `3be6fa8b03d544367b3b61de058de4784885aa20e7e26b8b95813e35dfe7e9f4` |
| `d10_4.cert` (uncompressed, 3,414,908 bytes) | `7dc38b37e9af02340fe0282b523452968e6d0e72320eff4c8840d9a67b4e7e6c` |
| `d10_4.vipr` (regenerated by `to_vipr.py`, 7,364,283 bytes; not included) | `f2e8a0650842e7e445a0fad2988b934aecf8731ee5ed3b200847f6e419173a8e` |

## 7. AI-use disclosure

This result was found with Claude (Anthropic, Claude Opus 5.5), directed by Dennis Koshta. The same AI system did all of
the following:
- formulated the problem;
- wrote the search, certificate and checking code;
- produced the certificate;
- ran the checks;
- did the literature search;
- wrote this note.

The `viprchk` checker is third-party software. The verifier `verifier/degree_sensitivity.py` was written by a separate
AI agent instance. **No human has checked the proof, the code or this note.**

## References

- [T37] T. Tao et al., *Optimization constants*, entry 37a "The degree–sensitivity exponent",
  github.com/teorth/optimizationproblems, `constants/37a.md` (commit 46dbf03, 2026-10-17).
- [NS94] N. Nisan, M. Szegedy, On the degree of Boolean functions as real polynomials, *Comput. Complexity* 4 (1994)
  301–313.
- [NW95] N. Nisan, A. Wigderson, On rank vs. communication complexity, *Combinatorica* 15 (1995) 557–565, footnote 1.
- [HKP11] P. Hatami, R. Kulkarni, D. Pankratov, Variations on the sensitivity conjecture, *Theory of Computing
  Library, Graduate Surveys* 4 (2011), Example 5.4.
- [Tal13] A. Tal, Properties and applications of Boolean function composition, *ITCS 2013*, 441–454.
- [P21] N. V. Proskurin, On separation between the degree of a Boolean function and the block sensitivity,
  arXiv:2101.08600 (2021), §5 and Theorem 11.
- [W22] J. Wellens, Relationships between the number of inputs and other complexity measures of Boolean functions,
  *Discrete Analysis* 2022:20; arXiv:2005.00566 (2020), Table 1 and Fact 3.11.
- [MMST] S. Maitra, C. S. Mukherjee, P. Stănică, D. Tang, On Boolean functions with low polynomial degree and higher
  order sensitivity, arXiv:2107.11205 (2021); *Discrete Appl. Math.* 365 (2025) 138–159.
- [Bal21] K. Balodis, Several separations based on a partial Boolean function, arXiv:2103.05593 (2021).
- [LZ26] Y. Li, J. Zhang, arXiv:2610.17006 and arXiv:2610.17100 (2026).
- [VIPR] K. K. H. Cheung, A. Gleixner, D. E. Steffy, Verifying integer programming results, *IPCO 2017*, LNCS 10328,
  148–160; checker at github.com/ambros-gleixner/VIPR.
