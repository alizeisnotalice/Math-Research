# upper004：谱二阶矩 uniform cube 下界审计

## Verdict

**接受，条件为已锁定的完整 fan/Fourier-density identity。** 该方法确实给出
整个 \(\|x\|_\infty\le1\) 上的严格 Green 下界

\[
G(x)\ge
[6.35027037455189259211051115037502266441401464984175997804459975185723667221633809207702364
\cdot10^{-8}\ \mathbin{\pm}\ 1.79\cdot10^{-98}].
\]

这是 uniform \(g\) 的下界，因而可写 \(g\ge\) 该区间下端；不是 \(I\) 的上界、
不是成本界，也不是 \(C_6\) 的新上界或下界。

## Exact reconstruction

Actual independent run (exit 0) used only the locked 21-cell exact fan and
python-flint at 512 bits:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
      /Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/upper004/independent_spectral_moment_audit.py

It independently obtains

\[
s_0=\frac{520066141}{2247297614217216000},\qquad
\sum_R\alpha_R={\rm fullMass}/46080,
\]

and

\[
E=\mathbb E_{\mu_1}\|\xi\|_2^2
=\frac{41159997228188655673}{685874959442940960000}.
\]

The full exponential mass is exactly

\[
46080s_0=\frac{520066141}{48769479475200},
\]

matching the fan's locked full_exp_minus_symbol_integral. Consequently

\[
C_0=\frac{46080s_0}{(2\pi)^6}
=1.73313043723422260282186019776173274879353387242806623199838969531188818014
\cdot10^{-10}.
\]

## Formula audit

On a base cone, write \(\xi=\sum_j u_jr_j\), where under the normalized
\(e^{-\psi}\) spectral measure the \(u_j\)'s are independent exponential
variables with rates \(q_j\). Therefore, for each coordinate \(i\),

\[
\mathbb E\xi_i^2=
\left(\sum_j{r_{ji}\over q_j}\right)^2+
\sum_j{r_{ji}^2\over q_j^2}.
\]

Summing this identity over \(i\) gives exactly the implemented expression

\[
\sum_i\left[\left(\sum_jr_{ji}/q_j\right)^2+
\sum_j(r_{ji}/q_j)^2\right].
\]

Multiplication by \(\alpha_R=|\det R|/\prod_jq_j\), followed by division by
\(s_0\), gives the base average. Signed permutations preserve Euclidean norm,
so the same quotient is the full \(W\)-average. More strongly, summing the
second-moment matrix over \(W\) kills every off-diagonal entry (an individual
sign flip changes its sign) and makes the six diagonals equal (permutations).
The normalized full covariance is consequently \(E I_6/6\).

Homogeneity gives a normalized \(t\)-spectral variable with second moment
\(t^{-2}E\). Because the full density is unconditional,

\[
\mathbb E[(x\cdot\xi)^2]
=\sum_i x_i^2\mathbb E\xi_i^2
\le E/t^2\quad(\|x\|_\infty\le1).
\]

Using the Fourier representation, its imaginary part cancels by evenness and
\(\cos u\ge1-u^2/2\), hence

\[
p_t(x)\ge C_0t^{-6}\left(1-\frac{E}{2t^2}\right).
\]

For \(T=\sqrt{E/2}\), the right side is nonnegative on \(t\ge T\). Since
\(p_t\ge0\), discarding \(0<t<T\) is allowed, and exact integration gives

\[
G(x)\ge C_0\int_T^\infty t^{-6}\left(1-\frac{E}{2t^2}\right)dt
=C_0\left(\frac1{5T^5}-\frac{E}{14T^7}\right).
\]

This proves the claimed entire-cube scope; it does not use the invalid
round003 floating direct-kernel samples.

## Limitations

The proof inherits only the previously audited full fan/Fourier-density
identity. It has no face integral estimate and its scale
\(6.35\cdot10^{-8}\) alone cannot improve a cost comparison. No \(C_6\)
breakthrough is asserted.

## Hashes

| artifact | SHA-256 |
|---|---|
| locked fan | 18d8aa578f84283d00e1c9c30522fb14cc49812092afd5716f2cc9ae68ba2798 |
| producer source | 872bc37259f5caac7da0146e825cb1d6ffa89439661ab63ee216a72434b07e43 |
| producer result | 159fa7e3cdd93b11568e3548008ff6d75be53759ecb6dc588cbd898cc7ee71e0 |
| independent checker | 005ae3c8aa0b4e5548f50d8093527b73836cac1a4a65d0546231a9ae5be3210e |
