#!/usr/bin/env python3
"""Check a branch-and-bound non-existence certificate (format in certify.py). Standard library only, exact integers.

Claim checked: there is no f:{0,1}^n -> {0,1} with deg f <= d (real multilinear degree), f(0) = 0 and f(e_i) = 1
for all i. The checker rebuilds the Mobius rows itself from (n, d); the certificate only supplies the branching
tree and integer multipliers.

    python3 -I check_cert.py FILE.cert
"""
import hashlib
import sys


def main(path):
    data = open(path, "rb").read()
    lines = data.decode("ascii").split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    head = lines[0].split()
    if head[:2] != ["p", "cert-degsens"] or len(head) != 4:
        raise SystemExit("bad header")
    n, d = int(head[2]), int(head[3])
    if not (2 <= n <= 16 and 1 <= d < n):
        raise SystemExit("n out of range")
    N = 1 << n
    pc = [bin(x).count("1") for x in range(N)]
    rows = [S for S in range(N) if pc[S] > d]
    # row r: coefficient (-1)^{|S|-|T|} on f_T for T subseteq S
    lo = [0] * N
    hi = [1] * N
    hi[0] = 0
    for i in range(n):
        lo[1 << i] = 1
    pos = 1
    leaves = 0
    maxdepth = 0
    # iterative preorder walk: stack of (mask, stage)
    stack = []

    def leaf(tokens):
        k = int(tokens[1])
        if len(tokens) != 2 + 2 * k:
            raise SystemExit("bad leaf")
        c = {}
        seen = set()
        for t in range(k):
            r, y = int(tokens[2 + 2 * t]), int(tokens[3 + 2 * t])
            if not (0 <= r < len(rows)) or r in seen:
                raise SystemExit("bad row index")
            seen.add(r)
            S = rows[r]
            T = S
            while True:
                c[T] = c.get(T, 0) + (y if (pc[S] - pc[T]) % 2 == 0 else -y)
                if T == 0:
                    break
                T = (T - 1) & S
        mx = 0
        for j, cj in c.items():
            mx += cj * (hi[j] if cj > 0 else lo[j])
        if not mx < 0:
            raise SystemExit("leaf certificate fails")

    while True:
        if pos >= len(lines):
            raise SystemExit("truncated")
        tok = lines[pos].split()
        pos += 1
        if tok[0] == "b":
            j = int(tok[1])
            if not (0 <= j < N) or lo[j] != 0 or hi[j] != 1:
                raise SystemExit("branch on fixed/invalid variable")
            hi[j] = 0
            stack.append(j)
            maxdepth = max(maxdepth, len(stack))
            continue
        if tok[0] != "l":
            raise SystemExit("bad line")
        leaf(tok)
        leaves += 1
        # pop finished 1-branches, then switch the innermost 0-branch to its 1-branch
        while stack and lo[stack[-1]] == 1:
            j = stack.pop()
            lo[j], hi[j] = 0, 1
        if not stack:
            break
        j = stack[-1]
        lo[j], hi[j] = 1, 1
    if pos != len(lines):
        raise SystemExit("trailing data")
    print(f"VALID: no f on {n} variables with deg <= {d} fully sensitive at 0 "
          f"({leaves} leaves, depth {maxdepth}, sha256 {hashlib.sha256(data).hexdigest()})")


if __name__ == "__main__":
    main(sys.argv[1])
