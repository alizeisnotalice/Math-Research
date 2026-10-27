# R034 q=1/20: improved explicit sufficient threshold

2026-10-28. **ACCEPT.** For every integer n with L=log n>=2^22, the accepted
R034 example q=1/20 implies

```
C_n<T-n/(40sqrt L),    T=nL-n log L.
```

Write l=log L. The accepted scalar expression is

```
W_n=T+(3l+1)n/L+(9l^2+1)n/L^2.
```

To prove W_n<T+n/(40sqrt L), it suffices to bound

```
H(L)=40[(3log L+1)/sqrt L
         +(9(log L)^2+1)/L^(3/2)]<1.
```

All four terms are decreasing in this range. The derivative signs for
log L/sqrt L and (log L)^2/L^(3/2) are determined by
1-(log L)/2<0 and 2log L-(3/2)(log L)^2<0 respectively; the pure powers
are also decreasing.

At L0=2^22, sqrt L0=2048 and log L0=22log2<77/5, since

```
exp(7/10)>1+7/10+(7/10)^2/2+(7/10)^3/6
          =12013/6000>2.
```

Thus the base value has the exact rational upper bound

```
H(L0)<40[(236/5)/2048+(53386/25)/2048^3]
      =2474666053/2684354560<1.
```

Monotonicity gives the inequality for every L>=L0. The R034 example applies
already at L>=4096, so it applies throughout this larger threshold range:

```
C_n<W_n-n/(20sqrt L)
   <T-n/(40sqrt L).
```

`python3 check_small.py` actually verified both displayed rational
comparisons using Fraction and wrote `small_receipt.json`. This is an
explicit threshold proof, not numerical searching or a statement that
2^22 is the smallest possible sufficient threshold. The R034 source hash
is recorded in that receipt. No frozen R033 file was modified.
