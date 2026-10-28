# R082 tau-tail explanation repair

The frozen producer README correctly states

    y S16'(y)-(5/2)S16(y)>0 for y>=4,

but its explanation that the cubic contribution alone dominates the
negative degrees 0,1,2 on this entire domain is false. The partial
polynomial through degree three at y=4 equals

    -5/2-3y/2-y^2/4+y^3/12 = -43/6 <0.

This does not invalidate the displayed S16 derivative conclusion.
Keep the positive quartic as well. All higher coefficients in
yS16'-(5/2)S16 are nonnegative, hence for y>=0

    yS16'-(5/2)S16
      >=p(y)=-5/2-3y/2-y^2/4+y^3/12+y^4/16.

With u=y-4>=0, exact expansion gives

    p(4+u)=53/6+(33/2)u+(27/4)u^2+(13/12)u^3+u^4/16>0.

Every coefficient is positive. This proves the required domain y>=4,
and therefore decrease of L^(5/2)/S16(L/50) wherever L/50>=4.
The registered base L0=2^24 satisfies this condition directly.

The independent exact-arithmetic check verifies the coefficients and
the countervalue -43/6 for the incomplete explanation. The producer
README, program, registration, receipt and their hashes remain untouched.
No new identity is assigned to that receipt. The R081 eventual proof
uses exponential limits and never used this erroneous cubic explanation;
its frozen hash and original eventual scope are unchanged.

The independent full R082 proof also uses the simpler valid inequality
exp(y)>=y^8/8! to bound tau. Its whole-tail envelope is a constant
times L^(-11/2), so acceptance does not rest on the repaired prose.
