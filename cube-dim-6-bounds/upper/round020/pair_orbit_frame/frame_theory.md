# A real body32 plus pair-orbit Cauchy frame

Let V be the32 unoriented cube-body representatives and
P={e_i+e_j,e_i-e_j:1<=i<j<=6}. Use the62 physical columns
A_beta=[V,beta P], beta>0, with an independent standard Cauchy variable
of common scale t on every column. This fixes one actual probability law,
not independently chosen coordinate times. The sole numeric weight in
this task is beta=1.

## Pair symbol, Gram, and unchanged canonical fan

For real a,b,
|a+b|+|a-b|=2max(|a|,|b|). Thus for z=sort_desc(|xi|),

    psi_pair(xi)=sum_(i<j)(|xi_i+xi_j|+|xi_i-xi_j|)
       =2sum_(i=1)^6(6-i)z_i.

On the whole closed sorted chamber xi1>=...>=xi6>=0, including all ties
and zero coordinates, this is the single linear functional
ell.xi, ell=(10,8,6,4,2,0). The new sum/difference hyperplanes lie on
chamber boundaries, not inside any old canonical cone. Therefore the old
21-cone body fan remains a valid geometric fan; only its rates change.

The true symbol and Gram are

    psi_beta=psi_body+beta psi_pair,
    A_beta A_beta^T=(32+10beta^2)I.

The latter follows because each support pair contributes
2(e_i e_i^T+e_j e_j^T), and every coordinate occurs in five pairs.
The spectral increment is linear in beta; the Gram increment is quadratic.
At beta=1 both the coordinate marginal rate and Gram scalar happen to42,
but their different beta dependencies cannot be interchanged.

For each existing integer ray r with old rate q=psi_body(r), the exact
new ray rate is

    q_new=q+beta*(10r1+8r2+6r3+4r4+2r5).

At beta=1 it is integral. Keep the ray vectors and their determinants
unchanged. Per cone the Fourier weight becomes det(rays)/product(q_new).
The old geometric coverage proof and exact ray-scaling bridge are frozen
inputs. No new arrangement or cone enumeration is necessary. Full signed
permutation symmetry is preserved by the pair orbit, so the old46080
group images remain the correct whole-space expansion.

The exported fan has a new hash and represents a different model. Old
kernel values, gradients, minima, source catalogs and overflow/genericity
certificates cannot be labeled as new-model certificates. The linear forms
x.r=0 use the same rays, but divided-difference collision forms involving
q_new can change. Subsequent evaluators must bind the new fan and recheck
their own exact input/range conditions.

## Positivity and the exact common-T/full-scale Markov interface

Independent convolution gives a symmetric probability law with Fourier
transform exp(-t psi_beta). It is a positive, mass-preserving, self-adjoint
convolution semigroup, contracting both L1 and L-infinity. The Gram identity
and l1/l2 comparison give

    sqrt(32+10beta^2)||xi||<=psi_beta(xi)
       <=(32+beta sqrt(220))||xi||.

The upper bound uses the old body32 L2 moment and the Euclidean norm of
ell. Hence the multiplier is integrable for every t>0 and its inverse is
a continuous density. It is strictly positive by selecting an invertible
six-column subset: its product Cauchy pushforward is strictly positive,
and convolution with the other probability factors preserves positivity.
The first-order scaling makes these kernels an L1 approximate identity,
so P_t is strongly continuous on L1. The L2 generator is-psi_beta(D),
with domain H1 and form domain H^(1/2). Test functions have Lphi in L1
and L-infinity by the finite sum of Cauchy second-difference integrals.

For K_T=integral_0^T p_t, positivity/Tonelli give mass exactly T, and
testing the semigroup derivative gives

    L_beta K_T=p_T-delta_0.

This is a distributional identity, not an assertion K_T lies in the L2
generator domain. For mass-preserving spatial dilation,

    (K_T)_r(x)=r^-6 K_T(x/r)=r^-1 K_(rT)(x),
    L_beta(K_T)_r=(p_(rT)-delta_0)/r.

The endpoint is rT. A different inward-point identity is
K_T(r y)=r^-5 K_(T/r)(y)>=K_T(y),0<r<=1. Both are for this same fixed
frame and one common radius r in all six coordinates.

The positive Markov average maximal inequality has weak-L1 constant1
for this semigroup. A short self-contained argument is as follows. For a
fixed step P=P_h, f>=0 in L1 and threshold a>0, define
R0=0 and Rn=(f-a+P R_(n-1))_+. These are increasing nonnegative L1
functions, Dn={Rn>0} has finite measure since
Dn subset{f+P R_(n-1)>a}, and integration on Dn gives
a|Dn|<=integral_Dn f. Positivity and P1=1 imply Rn dominates every
initial Cesaro sum minus its threshold. This proves the discrete average
weak bound1. Strong L1 continuity passes common-step Riemann averages
to finitely many positive rational times. Exhaust the rationals; Tonelli
and local absolute continuity of integral_0^s P_t f give one full-measure
set on which the average is continuous for every s>0. Thus all positive
times, and general L1 inputs via |f|, have the same bound1.

Consequently, if this new K_T is independently proved >=m>T/4096 a.e.
on the cube, the normalized kernel K_T/(64m) dominates the unit cube
average. Its dilation convolution is

    ((K_T/(64m))_r*f)=[T/(64m)]*(1/(rT))integral_0^(rT)P_t f dt.

The genuine six-dimensional centered common-radius maximal operator then
has weak constant at most T/(64m)<64 for Lebesgue measure, all positive
scales and general L1. No old m or singular barrier supplies the missing
new-model premise. This task certifies the interface, not its numerical
global lower-bound premise.

## Time normalization, with no time or point search

The coordinate marginal is exactly Cauchy with rate
s_beta=sum_columns|a_i,1|=32+10beta. A natural common time parameter is
tau=s_beta T. At beta=1 use T=tau/42. If main chooses a controlled shape
comparison at the old base T=6/125 marginal time tau=192/125, the new
physical cutoff is exactly T=32/875. This is one algebraic normalization,
not a predicted successful cutoff, a time grid, or an inherited point.

The old base T=.048 and old body-plus-axis(beta=1) T=2/45 have already
failed their own point tests. Their kernels are different from this pair
frame, so neither failure nor apparent Gram improvement certifies the
pair model. No kernel is evaluated at the proposed normalized time here.

## Band Taylor extension, not old-checker support

For a positive band B_(a,T), the coarea right inverse is
R=A_beta^T/(32+10beta^2), with ||Rh||2^2=||h||2^2/(32+10beta^2).
The Cauchy logarithmic second derivative bound-2/t^2 therefore gives

    D_h^2 B>=-2||h||2^2 J2_B/(32+10beta^2).

On a convex cell, Taylor's integral remainder has loss
D U2/(32+10beta^2). The original sparse six-body representation still
gives Phi(u)=exp(12asinh(u/12)): shift those body parameters and leave
the extra30 pair parameters fixed. This is a theoretical extension only.
Every band value/gradient, source record and singular lower barrier must
match the new frame. The frozen A015 checker permits only body/axes;
it has not been changed and does not presently approve pair-frame inputs.
