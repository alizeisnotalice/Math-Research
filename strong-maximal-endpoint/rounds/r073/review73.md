# Audit: exponential cap-two amplitude reduction

## Verdict

The proposed exponential-scale cap-two reduction is sound, conditional on
the archived high-amplitude splitting inequality that produces the three
displayed kernel terms in `exponential73.py/.wl`.  Its new constants create
a fixed exponential gap below the cap-two benchmark, and the core repair is
`1+O(Ne^{-N/100})`.  This is a **global supremum** reduction only; it does
not create the local fixed-reference mean hypothesis required by r072.

## Middle-range exponent

Let `b_4` be the positive root of

\[
 4(1-e^{-b_4/4})=1+e^{-b_4}.
\]

For cutoff `3/2`, set `a=\log(3/2)`, take `p=4` and
`s_N=e^{-N/100}`.  The lower middle-range logarithm in the checked formula
is

\[
 N(b_4+0.04)+3a-(N-1)\log(1+a).
\]

Its exponential rate is

\[
 b_4+0.04-\log(1+\log(3/2)).                             \tag{1}
\]

The rational inequalities recorded by the native verification suffice for
the required strict bound.  Namely, `e^{3/8}>16/11` gives
`e^{-3/8}<11/16`, hence the defining left side at `3/2` exceeds `5/4`,
while `e^{-3/2}<(11/16)^4<1/4`; therefore `b_4<3/2`.  Also
`e^{2/5}<3/2` and `e^{1/3}<7/5` give

\[
 \log(1+\log(3/2))>\log(7/5)>1/3.
\]

Thus

\[
 (1)<{3\over2}+{1\over25}-{1\over3}
 ={181\over150}<{5\over4}<\mu.                          \tag{2}
\]

The last inequality follows from `e^{5/4}<4`, since
`x\mapsto1+e^{-x}-x` is decreasing and is positive at `5/4`.
This establishes a fixed positive exponential-rate gap to the cap-two
scalar lower witness, whose logarithm is `N\mu-o(N)`.

## Other pieces and the repair

The upper middle endpoint is

\[
 N(b_4+0.04)+3\sqrt N-(N-1)\log(1+\sqrt N),
\]

which is `-\tfrac12N\log N+O(N)`.  The large/r=`1/2` kernel segment in the
checked formula has the same form after Stirling:

\[
 {N\over100}+N\log4-\log\Gamma(N)
 +(N-1)\{\log(\sqrt N+N(\log2+1/100))-
             \log(1+\sqrt N)\}
 =-\tfrac12N\log N+O(N).
\]

Both are superexponentially smaller than `T_{N,2}`.  Therefore the
high-amplitude coefficient is eventually below the repaired cap-two term.

For the retained core, with `r_N=2s_N`,

\[
 {3/2\over1-r_N}<2,\qquad
 (1-r_N)^{-N}=1+O(N e^{-N/100}).                           \tag{3}
\]

Consequently the usual subadditive split gives, for all sufficiently large
`N`,

\[
 {C_N\over T_{N,2}}\le(1-2e^{-N/100})^{-N}
 =1+O(N e^{-N/100}),                                      \tag{4}
\]

while the cap-two subclass gives the reverse inequality `C_N\ge T_{N,2}`.
So `C_N/T_{N,2}\to1` follows at the claimed rate.

## Scope / remaining gap

The arithmetic and exponent comparison are valid.  The numerical
Python/Wolfram output independently agrees with the constants but is only
a check; the proof still relies on the previously established analytic
high-amplitude splitting inequality.  Most importantly, (4) is a global
weak-ratio/supremum statement.  It does not say that an arbitrary source can
be clipped to cap two with `o(\log N/N)` loss at a particular retained
receiver.  The local reference/anchor obstruction in `reference73.md`
therefore remains untouched.


## Integration

The final amplitude73.md now displays the exact middle/extreme coefficients and disjoint-source maximum inequality explicitly, with r001 as their analytic source. The strict constant comparison also includes its rational polynomial step.
