**🇬🇧 English** · [🇫🇷 Français](README.fr.md)

# Goppa Codes and Lagrange Interpolation

> Master 2 thesis in Applied Mathematics (Numerical Computing), **University of Antananarivo**, 2026. Supervisor: Prof. Harinaivo ANDRIATAHINY.

Goppa codes are error-correcting codes built from finite fields and polynomials. They are the foundation of the **McEliece cryptosystem**, one of the main candidates for **post-quantum cryptography**. This project studies them and their link with **Lagrange interpolation**, and provides a small **from-scratch Python implementation** (no dependency) that you can run and test.

## The idea in 30 seconds

Let `L = {γ₁, …, γₙ}` be a set of distinct elements of a finite field and `g(x)` a polynomial with `g(γᵢ) ≠ 0`. The **Goppa code** Γ(L, g) is the set of words `c` such that

```
R_c(x) = Σ cᵢ / (x − γᵢ) ≡ 0  (mod g(x))
```

Lagrange interpolation links a word to a polynomial: `c = (c₁,…,cₙ) → F(x) = Σ cᵢ Lᵢ(x)` with `F(γᵢ) = cᵢ`. For a received word `r = c + e`, the **syndrome** `S(x) = Σ rᵢ / (x − γᵢ) mod g(x)` depends only on the error `e`, which is the starting point of decoding algorithms.

Main results presented in the thesis:
- dimension `k ≥ n − m·t` and minimum distance `d ≥ t + 1`, with `t = deg g`;
- for **binary irreducible** Goppa codes: `d ≥ 2t + 1` (parameters `[n, ≥ n − mt, ≥ 2t + 1]`);
- BCH codes appear as a special case (`g(z) = z^(d−1)`).

## What the code does

`src/goppa.py` implements, in pure Python:

| Part | Content |
|------|---------|
| Finite field | arithmetic in GF(2^m) (log/antilog tables) |
| Polynomials | addition, multiplication, evaluation, reduction modulo g |
| Goppa code | syndrome polynomial `S(x)`, parity-check matrix `H[i][j] = γⱼⁱ / g(γⱼ)` (expanded to bits), generator matrix (kernel of H), encoding |
| Decoding | syndrome table lookup, corrects up to `t` errors |
| Lagrange | interpolation over GF(2^m) and over the rationals |

Default example: **GF(8)** (from `x³ + x + 1`), `L` = the 8 field elements, `g(x) = x² + x + 1` (irreducible over GF(8)). This gives a binary code **[8, 2, 5]** that corrects **2 errors**.

> Scope: this is an educational implementation. The decoder uses a precomputed syndrome table. The efficient algebraic decoders (**Patterson**, **Berlekamp–Massey**, extended Euclid) are described in the thesis as future work.

## Run it

Requires Python 3.8+, no external package.

```bash
python src/goppa.py          # demo: encode, add 2 errors, decode, Lagrange
python tests/test_goppa.py   # run the tests (or: python -m pytest)
```

Example output (excerpt):

```
Binary Goppa code over GF(8): n=8, t=2, k=2, d=5
codeword c     : [1, 1, 1, 0, 1, 0, 1, 0]   S_c(x) = 0 (0 for a codeword)
received r     : [1, 0, 1, 0, 1, 0, 0, 0]   (errors at positions 1 and 6)
decoded        : [1, 1, 1, 0, 1, 0, 1, 0]   -> OK
F(L_i) == c_i for all i: True
```

## Tests

The test suite checks, among other things, that:
- GF(8) satisfies the field axioms;
- the binary parity-check matrix `H` defines **exactly the same code** as the definition `Σ cᵢ/(x − γᵢ) ≡ 0 (mod g)` (all 256 words);
- `H·Gᵀ = 0`, `k ≥ n − m·t` and `d ≥ 2t + 1`;
- **every** error pattern of weight ≤ 2 is corrected, for every codeword;
- Lagrange interpolation returns the right polynomial, e.g. the points (0, −1), (2, 3), (4, 15) give `x² − 1`.

## Thesis outline

1. **Mathematical foundations**: modular arithmetic, polynomials over Fq, finite fields, field extensions.
2. **Error-correcting codes**: linear codes, generator and parity-check matrices, dual code, cyclic codes, BCH codes.
3. **Goppa codes**: definition, parameters, parity-check matrix, binary Goppa codes, example over F₂³, asymptotic behaviour, decoding, use in cryptosystems.
4. **Goppa codes and Lagrange interpolation**: the interpolation problem, relation with Goppa codes, worked example, advantages, limits, perspectives, algorithm and Python implementation.

<!-- TODO: after the defense, add the link to the thesis PDF and slides here. -->

## Future work

- Full decoders: **Patterson** (binary), **Berlekamp–Massey**, **extended Euclid**.
- Faster finite-field operations.
- Cryptographic applications (McEliece / post-quantum cryptography).

## Repository structure

```
.
├── src/goppa.py            # implementation
├── tests/test_goppa.py     # tests
├── README.md               # English version
├── README.fr.md            # French version
└── LICENSE
```

## Author

**Haja Tahina RABARY** — University of Antananarivo, Faculty of Sciences, Mathematics and Computer Science.

## License

Released under the [MIT License](LICENSE).
