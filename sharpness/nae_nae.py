"""Write NAE_3 o NAE_3 (9 variables, degree 4, sensitivity 9) in the verifier's 2-line format. Standard library only.

    python3 -I sharpness/nae_nae.py > sharpness/nae_nae_9_4.txt

Coordinate i (1-based) of x is bit i-1 of x; blocks are (x1,x2,x3), (x4,x5,x6), (x7,x8,x9).
Bit x of the table is bit (3 - x mod 4) of hex digit x // 4.
"""


def nae(a, b, c):
    return 0 if a == b == c else 1


def f(x):
    bit = lambda i: (x >> i) & 1
    return nae(*(nae(bit(3 * j), bit(3 * j + 1), bit(3 * j + 2)) for j in range(3)))


n = 9
bits = [f(x) for x in range(1 << n)]
digits = [(bits[4 * k] << 3) | (bits[4 * k + 1] << 2) | (bits[4 * k + 2] << 1) | bits[4 * k + 3] for k in range(1 << (n - 2))]
print(f"{n} 4")
print("".join("0123456789abcdef"[d] for d in digits))
