# r019 independent theory audit — 2026-11-06

Accepted after read-only review of `rounds/r019/report.md`.

The exact split gives `B(q)=A(q+w)+F(q+w)-F(q)` and the displayed nonnegative
loss identity (8). Its derivative is `B'=-A(q+w)/n-f(q)<0`. The exact modular
factor is `K=1-exp(-w)+exp(-w)I_n`; for `w=L/2+log L`, this is
`1+sqrt(pi/2)/L+o(1/L)`, so `log K=o(1)`.

The r016/r017 tilted-density bounds make `f(q0+v)/f(q0)=1+o(1)` on the
required `O(L)` interval and give `A(q0+v)=1+O(L^2/n)`. Therefore the first
term of (8) is `w/n(1+o(1))`, the band term is `O(L^3/n^2)=o(L/n)`, and
`B(q0-r)-1=(r-1/2)L/n+o(L/n)` for fixed r. The resulting strict sign tests
prove `q_pl=q_soft-1/2+o(1)` without circularity.

The `q<q_pl` truncation/deep-layer realization is legitimate. Finally,
`a+(mu-1)/8=(9-mu)/24`; thus (12) has the correct log coefficient. No
blockers found. Existing 60/100-bit Wolfram data are numerical corroboration
only.
