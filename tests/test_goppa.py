import os, sys
from itertools import combinations, product
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from goppa import GF, GoppaCode, default_code, lagrange_gf, lagrange_rational, p_eval, p_mod, p_mul


def test_field_axioms():
    F = GF(3)
    for a in range(1, 8):
        assert F.mul(a, F.inv(a)) == 1
    for a, b, c in product(range(8), repeat=3):
        assert F.mul(a, F.add(b, c)) == F.add(F.mul(a, b), F.mul(a, c))


def test_g_has_no_root_in_support():
    code = default_code()
    assert all(p_eval(code.F, code.g, a) != 0 for a in code.L)


def test_inverse_of_x_minus_a_mod_g():
    code = default_code()
    for a, inv in zip(code.L, code.inv_cache):
        assert p_mod(code.F, p_mul(code.F, [a, 1], inv), code.g) == [1]


def test_parity_check_matches_goppa_definition():
    """Kernel of H (bits)  ==  { c : sum c_i/(x - L_i) = 0 mod g }  (all 256 words)."""
    code = default_code()
    for w in product([0, 1], repeat=code.n):
        w = list(w)
        in_kernel = all(s == 0 for s in code.syndrome(w))
        in_goppa = code.syndrome_poly(w) == []
        assert in_kernel == in_goppa


def test_generator_and_parameters():
    code = default_code()
    for row in code.G:
        assert all(s == 0 for s in code.syndrome(row))          # H . G^T = 0
    assert code.k >= code.n - code.F.m * code.t                  # k >= n - m t
    assert code.min_distance() >= 2 * code.t + 1                 # d >= 2t + 1


def test_corrects_all_patterns_up_to_t_errors():
    code = default_code()
    for c in code.codewords():
        for w in range(code.t + 1):
            for pos in combinations(range(code.n), w):
                r = c[:]
                for p in pos:
                    r[p] ^= 1
                assert code.decode(r) == c


def test_lagrange_gf_interpolates():
    code = default_code()
    for c in code.codewords():
        Fx = lagrange_gf(code.F, code.L, c)
        assert all(p_eval(code.F, Fx, a) == ci for a, ci in zip(code.L, c))


def test_lagrange_rational_thesis_example():
    coeffs = lagrange_rational([(0, -1), (2, 3), (4, 15)])
    assert [int(x) for x in coeffs] == [-1, 0, 1]                # x^2 - 1


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok  ", name)
