# R030 fine asymptotics: independent focused review

2026-10-27. **ACCEPT within the stated asymptotic scope.** For the already
defined scalar family U_n=inf_{t>1}F_n(t),

```
U_n=T+(3l+1)n/L+O(nl^2/L^2),
x_n=L-3l+(7l-2)/L+O(l^2/L^2)
```

hold as n tends to infinity, for every global minimizer x_n=n/t_n, where
L=log n, l=log L and T=nL-nl. All implicit constants can be chosen independent
of the minimizer. No uniqueness is needed. The new genuine-geometry operator
bound does not depend on these finer assertions.

## 1. Source and limits of review

The independently read source is
`/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round030/parameter_theory/proof.md`,
SHA256 `13efb4b0df7e6649ca631c14ccb17f3afcded6096061f0190cd45dfb60d12114`.
The focused scope is its Sections 1, 3 and 4: the Beta kernel, global
localization, uniform differentiation, fine optimal-parameter expansion, and
the resulting expansion of U_n. This is not certification of finite-dimensional
root completeness, a unique optimizer, or the separate endpoint root intervals.
Only this report was added; no source file was edited and no numerical search,
MATLAB, commit or push was performed.

We use the already proved original scalar family

```
A(x)=(1+x/n)^n,
B(x)=1-h^n+n(L-log x)-(n/2)I(x),
h=(1+x/n)/2,
I(x)=integral_0^H q_n(u)/(u+x/2) du,
H=(n-x)/2, q_n(u)=(1-u/n)^(n-1), 0<x<n,
F(x)=(sqrt(A(x))+sqrt(B(x)))^2.
```

The old coarse bound U_n=T+O(nl/L), and monotonic increase of B_n(t) in t,
are accepted dependencies. Only the upper part of that coarse asymptotic is
needed for the localization below; no fine optimizer location is assumed.

## 2. Beta identities and an explicitly uniform tail bound

Direct integration of powers gives, on the complete interval [0,n],

```
integral q_n=1,
integral u q_n=n/(n+1)<=1,
integral u^2 q_n=2n^2/[(n+1)(n+2)]<=2.
```

On [0,H] the mass and first moment are exactly

```
m=1-h^n,
r=n/(n+1)-n h^n+[n^2/(n+1)]h^(n+1).
```

These are unnormalized moments; no probability mass is lost silently. The
source's strict Jensen lower bound follows by normalizing q_n du by m, and
strictness follows from its positive density on the nondegenerate interval.

Writing a=x/2, the exact identity

```
1/(a+u)=1/a-u/a^2+u^2/[a^2(a+u)]
```

and the full moments imply

```
I(x)=2/x-4n/[(n+1)x^2]+R_n(x),
-2h^n/x<=R_n(x)<=16/x^3.                         (2.1)
```

For the lower remainder bound, the omitted tail integral is at most
h^n/a. For the upper remainder bound, the complete positive remainder is at
most a^(-3) integral u^2 q_n<=16/x^3. This proves both signs and constants.

For every x in [L/2,2L] and all sufficiently large n, 2L/n<=1/2 and hence
h<=3/4. This is a uniform statement over the full interval. More explicitly,
it already follows from n=e^L>=L^2/2>=4L when L>=18. The tails do not hide
constants depending on x: for every fixed integer M>=0,

```
h^n<=(3/4)^n<=exp(-n/4)<=M! 4^M n^(-M).           (2.2)
```

Here the final inequality follows from the M-th term of the positive
exponential series. The same bound, times at most 4/3, applies to h^(n-1).
Thus every fixed polynomial factor in n or L can be absorbed by choosing a
fixed M. This justifies all exponentially small-tail remainders below with
absolute constants and a common sufficiently-large-n threshold.

## 3. Moving endpoint: exact derivative and cancellation

For fixed n and interior x, ordinary Leibniz differentiation on the finite
smooth integration interval is legitimate. Since H'=-1/2 and
H+x/2=n/2,

```
I_x=-(1/2) integral_0^H q_n(u)/(u+x/2)^2 du
    -h^(n-1)/n.                                  (3.1)
```

The endpoint coefficient in the source is therefore correct. Importantly,
the endpoint term cancels exactly in the derivative of B:

```
-(h^n)_x=-(1/2)h^(n-1),
-(n/2)[-h^(n-1)/n]=+(1/2)h^(n-1),
B_x=-n/x+(n/4) integral_0^H q_n(u)/(u+x/2)^2 du.  (3.2)
```

This identity avoids even needing to discard a differentiated endpoint term
as asymptotically negligible. The elementary inequality

```
0<=a^(-2)-(a+u)^(-2)<=2u/a^3
```

and the full first moment give the explicit uniform error

```
|integral_0^H q_n(u)/(u+a)^2 du-4/x^2|
 <=4h^n/x^2+16/x^3.
```

Combining this with (2.2), for x in [L/2,2L],

```
B_x=-n/x+n/x^2+O(n/x^3).                         (3.3)
```

Similarly (2.1), including the bounded term 1-h^n, gives uniformly

```
B=n(L-log x)-n/x+O(n/x^2).                       (3.4)
```

The bounded term is absorbed because n/x^2 tends to infinity uniformly on
this interval. If needed independently, (3.1) and the same integral estimate
also give

```
|I_x+2/x^2|<=2h^n/x^2+8/x^3+h^(n-1)/n.
```

Thus the author's differentiation and uniform error estimates are sound.

## 4. Global localization is independent of the desired fine parameter

For sufficiently large n there is an interior global minimizer: the already
proved explicit candidate gives U_n<=T+O(nl/L)<2^n, whereas F approaches
2^n at t down to 1 and diverges at t tending to infinity. Continuity then
places a minimizer in a compact interior interval. No finite-n existence
claim is needed to establish the n-to-infinity result.

Let x_n be any such global minimizer. If x_n<=L/2, monotonicity of B in t and
the inherited general integral upper bound I<=2t/(n-1) imply

```
F(x_n)>=B( L/2 )>=T+n log 2-4n/L.
```

In this display B is written as a function of x=n/t. The positive n log 2
gap eventually exceeds the fixed coarse upper remainder K n l/L, ruling out
this region for every minimizer. If x_n>=2L, monotonicity of A gives

```
A(x_n)>=(1+2L/n)^n
       >=exp(2L-2L^2/n)>=n^2/e,
```

again contradicting the same coarse upper bound for large n. Hence all
global minimizers lie in [L/2,2L], where Sections 2-3 are already uniform.

From A(x_n)<=U_n<=2nL and
log A>=x-x^2/(2n), we obtain

```
x_n<=L+l+log 2+2L^2/n.
```

From B(x_n)>=n(L-log x_n)-2n/x_n and
B(x_n)<=U_n<=T+K n l/L, we get

```
log x_n>=l-Kl/L-2/x_n>=l-(Kl+4)/L,
x_n>=L exp[-(Kl+4)/L]=L-O(l).
```

Therefore x_n=L+O(l) for every global minimizer. The constants in these
inequalities are inherited only from a coarse bound and are independent of
the selected minimizer. This localization neither assumes nor uses the
later fine parameter expansion.

## 5. Unsquared stationarity and the coefficient 7l-2

Put alpha=A_x/A=n/(n+x)>0. At an interior minimizer,

```
F_x=(sqrt A+sqrt B)(A_x/sqrt A+B_x/sqrt B)=0,
sqrt(AB)=-B_x/alpha>0.                           (5.1)
```

B_x<0 follows either from the inherited monotonicity or directly from (3.2)
and x>1 in the localization region. Thus logarithms are legal, and the
stationarity identity is not obtained by a sign-ambiguous squaring step.
From (3.3), uniformly at every minimizer,

```
-B_x/alpha=(n/x)[1-1/x+O(x^(-2))](1+x/n).         (5.2)
```

Using only x=L+O(l), the logarithmic forms are

```
log A=x+O(L^2/n),
log B=L+l+O(l/L),
2log(-B_x/alpha)=2L-2l+O(l/L).
```

Substitution in log A+log B=2log(-B_x/alpha) first proves

```
x_n=L-3l+O(l/L).                                (5.3)
```

One may now legitimately bootstrap. Relation (5.3) gives

```
log x_n=l-3l/L+O(l^2/L^2),
log B=L+l-l/L+O(l^2/L^2),
2log(-B_x/alpha)=2L-2l+(6l-2)/L+O(l^2/L^2).
```

The last expression follows by taking the logarithm in (5.2): the term
-2log x supplies 6l/L and the factor 1-1/x supplies -2/L. Terms from
log(1+x/n) and log A-x are much smaller than l^2/L^2 by the elementary
exponential-versus-power bound. Consequently the exact logarithmic identity
yields

```
x_n=L-3l+(7l-2)/L+O(l^2/L^2).                   (5.4)
```

In particular, the additional +l/L comes from transferring the -l/L term
in log B to the other side; the reported coefficient 7 is correct.

## 6. Value at any minimizer

Equations (3.4) and (5.3) give

```
B=T+(3l-1)n/L+O(nl^2/L^2),
2sqrt(AB)=2n/L+O(nl/L^2),
A=O(n/L^3).
```

The second estimate can use the exact stationarity equation (5.1), rather
than estimating the square root independently. Summing these three terms
proves

```
U_n=T+(3l+1)n/L+O(nl^2/L^2).
```

All estimates were uniform on a common region containing every global
minimizer. Thus neither the value asymptotic nor (5.4) assumes a unique
minimizer. These are properties of the old scalar family U_n. They supply
no corresponding lower bound or asymptotic equality for the true operator
constant C_n, whose new geometric upper bound is a separate result.

## 7. Actual bounded symbolic checks

Using the existing local dependency directory named in the source's
reproduction command, SymPy 1.14.0 was actually run on this computer. It
returned `FINE_ASYMPTOTIC_EXACT_CHECKS_PASS` after verifying:

- the exact moving-endpoint coefficient -h^(n-1)/n;
- its exact cancellation with the derivative of -h^n in B_x;
- the derivative of the stated truncated first-moment expression;
- coefficient matching in the logarithmic stationarity identity, yielding
  7l-2.

No numerical optimizer, root search or MATLAB was run. These symbolic checks
support algebra and transcription; the uniform tail estimates, global
localization and asymptotic conclusions are proved above. There is no
remaining gap in the reviewed fine-asymptotic claim.
