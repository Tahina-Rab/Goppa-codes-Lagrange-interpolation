"""
Binary Goppa codes and Lagrange interpolation, from scratch (pure Python).

Educational implementation accompanying the Master 2 thesis
"Code de Goppa et interpolation de Lagrange" (University of Antananarivo).

Example used by default:
    field    : GF(2^3), built from the primitive polynomial x^3 + x + 1
    support  : L = all 8 elements of GF(8)
    Goppa g  : g(x) = x^2 + x + 1  (irreducible over GF(8), so g(a) != 0 for all a in L)
    code     : binary Goppa code Gamma(L, g)  with parameters [8, k, d], t = deg g = 2,
               d >= 2t + 1 = 5  ->  corrects up to 2 errors.

What is implemented
    * arithmetic in GF(2^m) and polynomials over GF(2^m)
    * Goppa syndrome polynomial   S(x) = sum_i r_i / (x - L_i)  mod g(x)
    * parity-check matrix         H[i][j] = L_j^i / g(L_j)  (expanded to bits)
    * generator matrix (kernel of H), encoding, syndrome decoding (table lookup)
    * Lagrange interpolation, over GF(2^m) and over the rationals

Not implemented: Patterson / Berlekamp-Massey decoders (listed as future work).
"""
from fractions import Fraction
from itertools import combinations, product


# ----------------------------------------------------------------------------
# 1. Finite field GF(2^m)
# ----------------------------------------------------------------------------
class GF:
    """GF(2^m) with elements stored as integers 0 .. 2^m - 1 (polynomial basis)."""

    def __init__(self, m=3, prim_poly=0b1011):  # x^3 + x + 1
        self.m, self.q, self.prim = m, 1 << m, prim_poly
        self.exp, self.log = [0] * (2 * self.q), [0] * self.q
        x = 1
        for i in range(self.q - 1):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & self.q:
                x ^= prim_poly
        for i in range(self.q - 1, 2 * self.q):
            self.exp[i] = self.exp[i - (self.q - 1)]

    def add(self, a, b):
        return a ^ b          # characteristic 2: addition = subtraction = XOR

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def inv(self, a):
        if a == 0:
            raise ZeroDivisionError("0 has no inverse")
        return self.exp[(self.q - 1) - self.log[a]]

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def pow(self, a, n):
        if n == 0:
            return 1
        if a == 0:
            return 0
        return self.exp[(self.log[a] * n) % (self.q - 1)]

    def bits(self, a):
        """Expand a field element into m bits (most significant first)."""
        return [(a >> (self.m - 1 - i)) & 1 for i in range(self.m)]


# ----------------------------------------------------------------------------
# 2. Polynomials over GF(2^m): list of coefficients, lowest degree first
# ----------------------------------------------------------------------------
def p_trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def p_add(F, a, b):
    n = max(len(a), len(b))
    a, b = a + [0] * (n - len(a)), b + [0] * (n - len(b))
    return p_trim([F.add(x, y) for x, y in zip(a, b)])


def p_scale(F, a, s):
    return p_trim([F.mul(c, s) for c in a])


def p_mul(F, a, b):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] ^= F.mul(x, y)
    return p_trim(out)


def p_eval(F, a, x):
    r = 0
    for c in reversed(a):          # Horner
        r = F.add(F.mul(r, x), c)
    return r


def p_mod(F, a, g):
    a, g = p_trim(a), p_trim(g)
    inv_lead = F.inv(g[-1])
    while len(a) >= len(g):
        coef = F.mul(a[-1], inv_lead)
        shift = len(a) - len(g)
        for i, c in enumerate(g):
            a[shift + i] ^= F.mul(coef, c)
        a = p_trim(a)
    return a


def inv_linear_mod_g(F, a, g):
    """Inverse of (x - a) modulo g(x), for g(a) != 0.

    g(x) - g(a) = (x - a) q(x)  =>  (x - a) * q(x) = -g(a)  (mod g)
    =>  1 / (x - a) = -q(x) / g(a)   (and -1 = 1 in characteristic 2).
    """
    ga = p_eval(F, g, a)
    if ga == 0:
        raise ValueError("g(a) must be non-zero for every a in the support")
    # synthetic division of g by (x - a): quotient q, highest degree first
    coeffs = list(reversed(g))
    q, carry = [], 0
    for c in coeffs[:-1]:
        carry = F.add(c, F.mul(carry, a)) if q else c
        q.append(carry)
    q = list(reversed(q))
    return p_scale(F, p_trim(q), F.inv(ga))


# ----------------------------------------------------------------------------
# 3. Lagrange interpolation
# ----------------------------------------------------------------------------
def lagrange_gf(F, xs, ys):
    """Interpolating polynomial F(x) = sum_i y_i L_i(x) over GF(2^m), with F(x_i) = y_i."""
    result = []
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        if yi == 0:
            continue
        num, den = [1], 1
        for j, xj in enumerate(xs):
            if j != i:
                num = p_mul(F, num, [xj, 1])        # (x - xj) = (x + xj)
                den = F.mul(den, F.add(xi, xj))
        result = p_add(F, result, p_scale(F, num, F.div(yi, den)))
    return result


def lagrange_rational(points):
    """Classical Lagrange interpolation over Q. Returns coefficients (lowest degree first)."""
    n = len(points)
    coeffs = [Fraction(0)] * n
    for i, (xi, yi) in enumerate(points):
        basis, den = [Fraction(1)], Fraction(1)
        for j, (xj, _) in enumerate(points):
            if j != i:
                basis = [Fraction(0)] + basis                     # multiply by x
                for k in range(len(basis) - 1):
                    basis[k] -= xj * basis[k + 1]                 # ... minus xj * old poly
                den *= xi - xj
        for k in range(len(basis)):
            coeffs[k] += yi * basis[k] / den
    return coeffs


# ----------------------------------------------------------------------------
# 4. Goppa code
# ----------------------------------------------------------------------------
class GoppaCode:
    def __init__(self, F, support, g):
        self.F, self.L, self.g = F, list(support), p_trim(g)
        self.n, self.t = len(self.L), len(self.g) - 1
        self.inv_cache = [inv_linear_mod_g(F, a, self.g) for a in self.L]
        self.H = self._parity_check_bits()
        self.G = self._generator_from_kernel()
        self.k = len(self.G)
        self._table = None

    # -- Goppa syndrome polynomial  S(x) = sum r_i / (x - L_i)  mod g(x)
    def syndrome_poly(self, word):
        s = []
        for r, inv in zip(word, self.inv_cache):
            if r:
                s = p_add(self.F, s, inv)
        return s

    # -- parity-check matrix: H[i][j] = L_j^i / g(L_j), i = 0..t-1, expanded to bits
    def _parity_check_bits(self):
        F, rows = self.F, []
        for i in range(self.t):
            row = [F.div(F.pow(a, i), p_eval(F, self.g, a)) for a in self.L]
            for b in range(F.m):
                rows.append([F.bits(x)[b] for x in row])
        return rows

    def syndrome(self, word):
        return tuple(sum(h * w for h, w in zip(row, word)) % 2 for row in self.H)

    def _generator_from_kernel(self):
        """Basis of ker(H) over GF(2) (Gauss-Jordan elimination)."""
        M = [row[:] for row in self.H]
        n, pivots, r = self.n, [], 0
        for c in range(n):
            piv = next((i for i in range(r, len(M)) if M[i][c]), None)
            if piv is None:
                continue
            M[r], M[piv] = M[piv], M[r]
            for i in range(len(M)):
                if i != r and M[i][c]:
                    M[i] = [x ^ y for x, y in zip(M[i], M[r])]
            pivots.append(c)
            r += 1
        free = [c for c in range(n) if c not in pivots]
        basis = []
        for f in free:
            v = [0] * n
            v[f] = 1
            for row, pc in zip(M, pivots):
                if row[f]:
                    v[pc] = 1
            basis.append(v)
        return basis

    def encode(self, message):
        assert len(message) == self.k
        return [sum(m * g[i] for m, g in zip(message, self.G)) % 2 for i in range(self.n)]

    def codewords(self):
        return [self.encode(list(m)) for m in product([0, 1], repeat=self.k)]

    def min_distance(self):
        return min(sum(c) for c in self.codewords() if any(c))

    # -- decoding: syndrome table for all error patterns of weight <= t
    def _build_table(self):
        self._table = {}
        for w in range(self.t + 1):
            for pos in combinations(range(self.n), w):
                e = [1 if i in pos else 0 for i in range(self.n)]
                self._table[self.syndrome(e)] = e

    def decode(self, received):
        if self._table is None:
            self._build_table()
        e = self._table.get(self.syndrome(received))
        if e is None:
            return None                       # more than t errors: not correctable
        return [(r + x) % 2 for r, x in zip(received, e)]


def default_code():
    F = GF(3)
    return GoppaCode(F, support=list(range(F.q)), g=[1, 1, 1])     # g(x) = x^2 + x + 1


# ----------------------------------------------------------------------------
# Demo
# ----------------------------------------------------------------------------
def fmt_poly(F, p):
    if not p:
        return "0"
    names = {v: ("a^%d" % F.log[v] if v != 1 else "1") for v in range(1, F.q)}
    terms = []
    for i, c in enumerate(p):
        if c:
            terms.append((names[c] if (c != 1 or i == 0) else "") + ("x^%d" % i if i > 1 else "x" if i == 1 else ""))
    return " + ".join(reversed(terms))


def main():
    code = default_code()
    F = code.F
    print("Binary Goppa code over GF(8): n=%d, t=%d, k=%d, d=%d" % (code.n, code.t, code.k, code.min_distance()))
    print("Parity-check matrix H (bits):")
    for row in code.H:
        print("  ", row)
    print("Generator matrix G:")
    for row in code.G:
        print("  ", row)

    message = [1, 0] if code.k == 2 else [1] + [0] * (code.k - 1)
    c = code.encode(message)
    print("\nmessage        :", message)
    print("codeword c     :", c, " S_c(x) =", fmt_poly(F, code.syndrome_poly(c)), "(0 for a codeword)")

    r = c[:]
    for pos in (1, 6):                                  # two errors
        r[pos] ^= 1
    print("received r     :", r, "(errors at positions 1 and 6)")
    print("syndrome S(x)  :", fmt_poly(F, code.syndrome_poly(r)))
    print("decoded        :", code.decode(r), "-> OK" if code.decode(r) == c else "-> FAIL")

    Fx = lagrange_gf(F, code.L, c)
    print("\nLagrange polynomial F(x) of the codeword over GF(8):", fmt_poly(F, Fx))
    print("F(L_i) == c_i for all i:", all(p_eval(F, Fx, a) == ci for a, ci in zip(code.L, c)))

    pts = [(0, -1), (2, 3), (4, 15)]
    print("\nClassical Lagrange over Q for", pts, "->", [str(x) for x in lagrange_rational(pts)], "(= -1 + 0x + 1x^2)")


if __name__ == "__main__":
    main()
