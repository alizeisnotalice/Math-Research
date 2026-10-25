# Exact local partial fractions, real folding, and the Green integral

All formulas below are rational identities in an indeterminate T or z, not
evaluations at a positive time. Fix a six-node multiset B of real rational
nodes; zero and repetitions are allowed. Write its distinct nodes as b with
multiplicities m_b, summing to6, and set

`F_B(T)=prod_b(T-i*b)^(-m_b)`.

## 1. Local complex partial coefficients

At the pole T=i*b let s=T-i*b. For every other node c,

`T-i*c=i*(b-c+s/i)`.

Define rational coefficients a_k by the truncated Taylor series

`prod_(c!=b)(b-c+u)^(-m_c)=sum_(k=0)^(m_b-1) a_k*u^k+O(u^m_b)`.

The elementary negative-binomial factor is

`(d+u)^(-n)=sum_k (-1)^k binom(n+k-1,k) d^(-n-k)u^k`.

All d=b-c are nonzero by DISTINCT-node grouping; no zero-node assumption on
b itself is made. Substitution gives the principal part at i*b:

`sum_(k=0)^(m_b-1) i^(m_b-k-6)*a_k/(T-i*b)^(m_b-k)`.

With j=m_b-k, call a_k the real coefficient A_(b,j) before weighting.
Subtracting all these principal parts from F_B leaves a rational function
without finite poles. It is a polynomial. Since F_B and every principal term
vanish at infinity, that polynomial is zero. This proves the full partial
identity, including zero poles, repeated poles and sixfold poles.

The scalar checker verifies the identity independently by multiplying by
the degree6 complex denominator and obtaining exactly1, for the eight
registered cases. It does not sample T.

## 2. Whole-kernel aggregation and normalization

For the frozen C022 tuple rows, let w_B be the stored complete weight. Sum
`A_(b,j)=sum_B w_B*a_(B,b,j)` in exact fmpq arithmetic. C022 already included
stabilizer12 and both global signs; neither is added here. The identity is

`p_T(x0)=(2*pi)^(-6) Re sum_(b,j) i^(j-6) A_(b,j)/(T-i*b)^j`.

The global-negation merge preserves the real part of the density, not
necessarily its chosen complex representative. Accordingly do not assume a
signed-pole table is automatically conjugate symmetric.

For b>0, conjugating the negative-b denominator inside a REAL part gives

`C_(b,j)=A_(b,j)+(-1)^j*A_(-b,j)`.

This parity factor is necessary. Counting nonzero signed A entries before
this fold would not give the reduced real denominator.

At b=0, an odd j has a purely imaginary coefficient i^(j-6) and T^(-j) is
real: its real part is identically0. This algebraic fact is recorded; the
input terms were still accumulated. At even j=2,4,6 the real contributions
are respectively `A_(0,2)T^-2`, `-A_(0,4)T^-4`, `A_(0,6)T^-6`. Their distinct
orders must vanish separately for a bounded short-time density. B023 checks
each coefficient is exactly0, rather than imposing that value or dropping
zero-node tuples. C022's off-S3 regularity is consistent with the resulting
cancellation but is not substituted for the computation.

## 3. Real local numerator and the true denominator degree

For b>0 and z=T^2, expansion of (T+i*b)^j gives

`Re[i^(j-6)/(T-i*b)^j] = N_(b,j)(z)/(z+b^2)^j`,

where the ascending-power polynomial is

`N_(b,j)(z)=sum_(r=0)^floor(j/2) (-1)^(j-r-3)*binom(j,2r)*b^(j-2r)*z^r`.

Its degree is at most3, and strictly less than j. Evaluating at the pole,

`N_(b,j)(-b^2)=(-1)^(j+1)*2^(j-1)*b^j != 0`.

Indeed the signs of every even-binomial summand agree there and
`sum_r binom(j,2r)=2^(j-1)` for j>=1. For a fixed b, let J be the largest
order with C_(b,J)!=0. Its pole of order J at z=-b^2 cannot be cancelled by
lower orders; its leading local coefficient is C_(b,J)N_(b,J)(-b^2)!=0.
Other positive b values give distinct b^2, whose linear factors are coprime.
The true reduced denominator, up to a nonzero constant, is therefore

`prod_(b>0 with some C!=0)(z+b^2)^(max{j:C_(b,j)!=0})`.

Its degree is the sum of these local orders. This conclusion requires no
global numerator construction, polynomial gcd, or root isolation. The
optional ordinary principal table in z is obtained just by the local shift
`z=w-b^2` in N; each shift still has degree at most3. Independent saved-data
checks compare both representations and each highest-pole witness.

## 4. Exact integration: only complex-pole order1 contributes to G

All zero REAL poles are cancelled before integrating. Odd zero terms are
identically0. For b!=0 and j>=2,

`integral_0^infinity (T-i*b)^(-j)dT = (-i*b)^(1-j)/(j-1)`.

Multiplying by `i^(j-6) A_(b,j)` leaves a purely imaginary value: its power
of i is i^(-5) times a real sign. Hence the real integral is exactly0.
This integral is legitimate for j>=2: the complex term is bounded at0 and
integrable at infinity. No divergent complex integral at j=1 is used.

For j=1 directly take the real part:

`Re[-i*A_(b,1)/(T-i*b)] = A_(b,1)*b/(T^2+b^2)`.

Its real integral is `(pi/2)*sign(b)*A_(b,1)`. Every resulting nonzero-pole
real term is absolutely integrable on the positive half-line. The sum is
finite, so it may be integrated term by term after the exact zero-pole
cancellation. Thus with the correct dT measure,

`Gamma_G=sum_(b!=0) sign(b)A_(b,1)=sum_(b>0) C_(b,1)`,

`G(x0)=(2*pi)^(-6)*(pi/2)*Gamma_G=Gamma_G/(128*pi^5)`.

As an independent representation check, for the ordinary z-principal table
use

`integral_0^infinity (T^2+b^2)^(-j)dT = (pi/2)*binom(2j-2,j-1)/(4^(j-1)*b^(2j-1))`.

This follows from the j=1 arctangent integral and integration by parts,
`I_j=(2j-3)I_(j-1)/(2(j-1)b^2)`. The checker sums the rational coefficients
of this expression and compares them exactly to Gamma_G from the signed and
folded tables. No positive-time kernel is evaluated.

The symbolic Gamma_G is new complete kernel evidence, even though no T was
sampled. It is explicitly counted as part of ONE full symbolic computation.
It supplies no global cube lower bound or new C6 inequality. This first
stage ends after the partial tables, degree, Gamma_G and saved-data checks.
