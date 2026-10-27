# Proposal received for independent R030 audit

Provenance: parent-agent task message received on 2026-10-27; this snapshot is
the precise claim audited, not a separately discovered author manuscript.
Baseline clone HEAD: f65a5bf92911ebed4cf803071bd7693978e5e582.

For integer n with L=log(n)>=18, put l=log(L), T=nL-nl and
U_n=inf_{t>1}(sqrt(A_n(t))+sqrt(B_n(t)))^2, with
A_n(t)=(1+1/t)^n,
B_n(t)=1-h^n+n log(t)-(n/2)I,
h=(1+1/t)/2,
I=integral_1^t [1+1/(2t)-1/(2rho)]^(n-1) d rho/rho.

Proposed conclusion:

U_n <= T+(3l+1)n/L+(9l^2+1)n/L^2.

It strictly improves the previous upper expression
T+(6l+5)n/L on the previous common range L>=32.
The parameter is x=n/t=L-3l. The new integral estimate is
I>=2(1-2/n)/(x+2), obtained from the exact substitution
I=integral_0^((n-x)/2) (1-u/n)^(n-1)/(u+x/2) du,
the mass m=1-h^n, the first moment r<=1, and Cauchy-Schwarz.

No reduction of the leading coefficient 1, no lower bound for the true C_n,
and no global parameter optimality is claimed.
