# R083 fixed eventual-candidate constants

This is a single registered Fraction calculation of fixed constants. It does not establish a dimension threshold or a C_n theorem, and does not prove the local s=9/2 profile, target typicality, transport, weighted source-cap lemma, actual-history bridge or operator transfer. In particular, the choices below are arithmetic margins available to such proofs, not proof that those bounds hold at any finite dimension.

At s=9/2, a=(s+1)/(2s)=11/18 and mu=7/2+(11/2)log(11/18). Exact rational 64-term atanh-series enclosures give mu<0.7914<0.79143. If an independent argument establishes the 0.79143 x target event and its transport error at most 0.00001 x, then the transported 0.79144 x divided by (1-e) is 0.79144791448... x <0.79147 x. The 0.6 pin upper bound exceeds the log(18/11) enclosure, and chi-0.79147=0.00003. These arithmetic facts do not themselves establish the target event or derivative transport.

For qmax=9999/10000, trim=1/1000000, rlo=1+trim, rhi=s*qmax-trim=4.499549, and e=1/100000, the exact lower log enclosure verifies (1-e)log(rhi/rlo)>1.50396112>mL=1.503. With Ret=.99999 and the complete R081 P0, Ulo, Uhi, M and b formulas, (1-e)Gamma=1.50280462097... > gamma=1.501. The fixed weighted-cap parameters nu=qmin-(1+qmin)e and eta=(1-qmax)/(2qmax) are strictly positive; their interface still needs proof for this geometry.

The tau exponent is sigma/2-chi(1-sigma)=0.0001895>0. The log enclosure proves log(1000/387)<0.95=price; gamma-price-c=0.001 exactly for c=0.55. The complete scalar remainder 1/n+exp(-sigma L)+2 sqrt(L) exp(-sigma L/2) tends to zero because sigma=613/1000>0, so the fixed positive .001 allowance is eventually sufficient **if** the assumed scalar estimates apply. No finite payment or threshold was evaluated.

All 18 gates passed in actual PID 29125. Exact fractions, including log enclosures and margins, are in `result.json`. Registration preceded the run; no parameter grid, threshold run, shared MATLAB, commit or push was used. Reproduce with `python3 verify.py` (overwrites the PID-bearing result).
