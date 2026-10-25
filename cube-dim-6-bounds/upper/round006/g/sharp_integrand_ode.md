# Sharper valid coordinate ODE, and a withdrawn shortcut

Let N=32 and let A have the 32 body direction representatives as columns.
For one coordinate, c=A^T e_i/N has entries sigma_v/N with sigma_v in {+1,-1}.
At any coarea parameter point set

    F(s) = product_v f_t(z_v + sigma_v s/N),
    w_v = sigma_v z_v/t,
    a_v = -2w_v/(1+w_v^2),
    b_v = (1-w_v^2)/(1+w_v^2).

Then a_v^2=1-b_v^2, and the normalized second logarithmic derivative
of each Cauchy factor is -b_v(1+b_v). Consequently

    N^2 t^2 F''(0)/F(0)
      = (sum a_v)^2 - sum (b_v+b_v^2)
      <= N sum a_v^2 - sum (b_v+b_v^2)
      = N^2 - (N+1)sum b_v^2 - sum b_v
      <= N^2 + N/[4(N+1)].

The last step completes the square separately for every b_v. Hence the
coefficient is

    c_N = 1 + 1/[4N(N+1)],  c_32 = 4225/4224.

This is sharp for the unconstrained parameter integrand: choose all b_v
equal to -1/[2(N+1)] and all a_v with the same sign. Sharpness of this
parameter bound does not assert sharpness after the coarea integral.

The first two derivatives are dominated locally by a constant times F:
the original real multiplicative comparison handles nearby s, the score
is bounded by 1/t per factor, and each second logarithmic derivative has
absolute value at most 2/t^2. Differentiation under the coarea integral is
therefore legitimate. Positive integration gives

    partial_ii p_t <= (4225/4224)t^-2 p_t,
    partial_ii H_e <= (4225/4224)e^-2 H_e.

The already proved one-dimensional maximum principle then gives the same
positive tensor-product endpoint interpolant as round004, with

    k_e = sqrt(4225/4224)/e = 65/(8 sqrt(66)e),
    H_e >= L product_i cosh(k_e r_i)^-1

when every box vertex has lower bound L. This replaces the former
129/128 coefficient. At fixed endpoint loss, the allowed equal side lengths
improve by sqrt(4257)/65 = 1.003779838... . The improvement is strict but
small; it cannot close the full-face coverage problem with the current
finite budget. The actual eight local boxes use the independently audited
multiplicative inequality, not this new ODE estimate.

An initially considered coefficient `(33/512)(33/32)^32` is **withdrawn**.
The pure imaginary-axis tube estimate does not provide the required bound
on an entire Cauchy circle relative to p_t(x): its points have varying real
parts, and the density factor becomes p_t(x+Re(z)e_i). Thus applying the
Cauchy derivative estimate with the center density would be invalid. No
point or box acceptance depends on that shortcut.
