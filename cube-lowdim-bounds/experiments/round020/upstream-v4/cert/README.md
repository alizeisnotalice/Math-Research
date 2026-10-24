# Certificate format

A certificate is a JSON object with the following fields.

| field | type | meaning |
|---|---|---|
| `alpha` | rational string `"p/q"` | the exponent α of the operator, 1 < α < 2 |
| `R` | rational string | the truncation radius of the singular part |
| `h` | rational string | the spacing of the B-splines |
| `N` | integer | the coefficients a_ij vanish for max(\|i\|, \|j\|) > N |
| `c` | rational string | the coefficient of the singular part |
| `a` | object | maps `"I,J"` with 0 ≤ J ≤ I ≤ N to a rational string; missing keys mean 0 |

With r = |u| + |v|, the kernel is

    K(u, v) = c Ψ(r) + Σ_{i,j} a_ij B(u/h - i) B(v/h - j),        Ψ(r) = max(r^(-α) - R^(-α), 0),

where a_ij = a["max(|i|,|j|),min(|i|,|j|)"] for -N ≤ i, j ≤ N, and B is the centred cubic B-spline:
B(x) = (4 - 6x² + 3|x|³)/6 for |x| ≤ 1, B(x) = (2 - |x|)³/6 for 1 ≤ |x| ≤ 2, and B(x) = 0 otherwise.
So K is even and has the symmetries of the square.

For a bound C given on the command line, the checkers verify

- (C1) K ≥ 0 on the plane minus the origin;
- (C2) K ≥ 1 on the diamond |u| + |v| ≤ 1 minus the origin;
- (C3) (1/2) ∫K ≤ C, where ∫K = c · 2α R^(2-α)/(2-α) + h² Σ_{i,j} a_ij;
- (C4) AK(x) ≥ 0 for almost every x ≠ 0, where

      Af(x) = Σ_{j=1,2} ∫ [f(x + t e_j) + f(x - t e_j) - 2 f(x)] |t|^(-1-α) dt.

By Theorem 2.2 and Lemma 3.1 of the paper, (C1)-(C4) imply c_2 ≤ C. By the symmetries of K and of A it suffices
to check (C1), (C2) and (C4) on 0 ≤ v ≤ u.

`SHA256SUMS` lists the SHA-256 digests of the certificate files.
