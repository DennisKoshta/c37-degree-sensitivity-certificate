"""Branch-and-bound non-existence certificate for 'f:{0,1}^n -> {0,1}, deg f <= d, f(0)=0, f(e_i)=1'.

The system (checked by check_cert.py, stdlib only):
  variables f_T for every T subseteq [n] (mask 0..2^n-1), 0 <= f_T <= 1, f_0 = 0, f_{e_i} = 1;
  rows r = 0, 1, ...: the masks S with |S| > d in increasing order, sum_{T subseteq S} (-1)^{|S|-|T|} f_T = 0.
Certificate (text, preorder):
  "p cert-degsens n d"
  "b M"            branch on f_M: the f_M = 0 subtree follows, then the f_M = 1 subtree
  "l k r1 y1 ... rk yk"   leaf: integers y on rows; c = sum_r y_r A_r must satisfy max_{f in box} c.f < 0
                          (box = [0,1]^N with the root fixings and the branch fixings on the path), so A f = 0 is
                          infeasible there.
Every 0/1 point lies in exactly one leaf box, so all leaves infeasible => no such f.

    python -m problems.degree_sensitivity.certify N D OUT.cert [--rule mostfrac]
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np


def full_system(n: int, d: int):
    from scipy.sparse import coo_matrix
    N = 1 << n
    pc = np.array([bin(x).count("1") for x in range(N)])
    rows = [S for S in range(N) if pc[S] > d]
    r, c, v = [], [], []
    for i, S in enumerate(rows):
        T = S
        while True:
            r.append(i)
            c.append(T)
            v.append((-1) ** (pc[S] - pc[T]))
            if T == 0:
                break
            T = (T - 1) & S
    return rows, coo_matrix((v, (r, c)), shape=(len(rows), N)).tocsc()


class Certifier:
    def __init__(self, n, d, rule="mostfrac"):
        import highspy
        self.hs = highspy
        self.n, self.d, self.rule = n, d, rule
        self.rows, self.A = full_system(n, d)
        self.Ar = self.A.tocsr()
        N = 1 << n
        self.N = N
        self.lo, self.hi = np.zeros(N), np.ones(N)
        self.hi[0] = 0
        for i in range(n):
            self.lo[1 << i] = 1
        self.h = self._lp(self.A, slack=False)
        self.out = []
        self.nodes = self.leaves = 0

    def _lp(self, A, slack):
        hs = self.hs
        m, N = A.shape
        if slack:  # phase 1: A f + s+ - s- = 0, minimise sum s
            from scipy.sparse import hstack, identity
            I = identity(m, format="csc")
            A = hstack([A, I, -I]).tocsc()
        lp = hs.HighsLp()
        nc = A.shape[1]
        lp.num_col_, lp.num_row_ = nc, m
        cost = np.zeros(nc)
        if slack:
            cost[N:] = 1
        lp.col_cost_ = cost
        lp.col_lower_ = np.zeros(nc)
        up = np.full(nc, np.inf)
        up[:N] = 1
        lp.col_upper_ = up
        lp.row_lower_ = lp.row_upper_ = np.zeros(m)
        lp.a_matrix_.format_ = hs.MatrixFormat.kColwise
        lp.a_matrix_.start_, lp.a_matrix_.index_ = A.indptr, A.indices
        lp.a_matrix_.value_ = A.data.astype(float)
        h = hs.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("threads", 1)
        h.passModel(lp)
        return h

    def solve(self):
        idx = np.arange(self.N, dtype=np.int32)
        self.h.changeColsBounds(self.N, idx, self.lo, self.hi)
        self.h.run()
        if self.h.getModelStatus() == self.hs.HighsModelStatus.kInfeasible:
            return None
        return np.array(self.h.getSolution().col_value)

    def farkas(self):
        """Integer y with max over the box of (y^T A) f < 0, from the phase-1 LP's row duals."""
        if not hasattr(self, "h1"):
            self.h1 = self._lp(self.A, slack=True)
        idx = np.arange(self.N, dtype=np.int32)
        self.h1.changeColsBounds(self.N, idx, self.lo, self.hi)
        self.h1.run()
        y = np.array(self.h1.getSolution().row_dual)
        for sign in (1, -1):
            for K in (1 << 10, 1 << 16, 1 << 24, 1 << 32):
                yi = np.round(sign * y * K)
                nz = np.nonzero(yi)[0]
                yv = np.array([int(v) for v in yi[nz]], dtype=object)
                if self._check(nz, yv):
                    return nz, yv
        raise RuntimeError("could not round a Farkas certificate")

    def _check(self, nz, yv):
        c = [0] * self.N
        for r, yr in zip(nz.tolist(), yv.tolist()):
            row = self.Ar.getrow(r)
            for j, a in zip(row.indices.tolist(), row.data.tolist()):
                c[j] += int(a) * int(yr)
        mx = 0
        for j in range(self.N):
            if c[j] > 0:
                mx += c[j] * int(self.hi[j])
            elif c[j] < 0:
                mx += c[j] * int(self.lo[j])
        return mx < 0

    def run(self, x=None):
        self.nodes += 1
        x = self.solve()
        if x is None:
            nz, yv = self.farkas()
            self.leaves += 1
            self.out.append("l %d %s" % (len(nz), " ".join(f"{r} {y}" for r, y in zip(nz.tolist(), yv.tolist()))))
            return False
        frac = np.abs(x - np.round(x))
        free = (self.lo != self.hi)
        cand = np.where((frac > 1e-6) & free)[0]
        if cand.size == 0:
            raise RuntimeError("integral LP point: a solution exists?")
        j = int(cand[np.argmax(frac[cand])])
        self.out.append(f"b {j}")
        for v in (0, 1):
            olo, ohi = self.lo[j], self.hi[j]
            self.lo[j] = self.hi[j] = v
            self.run()
            self.lo[j], self.hi[j] = olo, ohi
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("d", type=int)
    ap.add_argument("out")
    a = ap.parse_args()
    sys.setrecursionlimit(100000)
    t = time.perf_counter()
    c = Certifier(a.n, a.d)
    c.run()
    with open(a.out, "w") as fh:
        fh.write(f"p cert-degsens {a.n} {a.d}\n")
        fh.write("\n".join(c.out) + "\n")
    print(f"nodes {c.nodes} leaves {c.leaves} sec {time.perf_counter() - t:.1f}")


if __name__ == "__main__":
    main()
