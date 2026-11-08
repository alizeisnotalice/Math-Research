# Exact radial envelope and global centroid coordinates

Fix b>=2 and c>0. All logarithms are natural. Let Z_i=E_i+cJ_i, where E_i are independent Exp(1), J_i are independent Bernoulli((1+e^c)^(-1)), and retain the complete joint label J. Put S=sum Z_i, H_b={y:sum y_i=0}, and Y=Z-S/b. Use Euclidean Hausdorff measure dH on H_b and counting measure on J, unless stated otherwise.

Define

\[
 Q(s)=\sum_{k=0}^b {b\choose k}(s-ck)_+^{b-1},\qquad
 H(s)=Q(s)/(b-1)!.
\]

The normalized conditional density is

\[
 k_s(y,J)=\frac1{\sqrt b H(s)}\mathbf1_{D_{s,J}}(y),
 \quad D_{s,J}=\{y:y_i\ge cJ_i-s/b\ \forall i\}.
 \tag{1}
\]

Indeed, the joint unconditioned density of (Z,J) is e^(-sum Z)/(1+e^(-c))^b on Z_i>=cJ_i: the Bernoulli factor e^(-c|J|) exactly cancels the shift in the exponential density. Coarea contributes 1/sqrt(b), and

\[
 f_{b,c}(s)=e^{-s}H(s)/(1+e^{-c})^b,
 \quad |D_{s,J}|_H=\sqrt b\,(s-c|J|)_+^{b-1}/(b-1)!.
\]

Uniform Haar probability on labels multiplies k_s by 2^b. Coordinate measure dy_1...dy_(b-1) removes the factor sqrt(b). No conclusions below depend on this choice if both kernel and measure are converted consistently.

## Exact envelope, including crossing knots

For 0<a<d<infinity the sets D_(s,J) are nested as s increases, while H(s) is strictly increasing and positive. A point already present at s=a attains its maximal density there. A point entering at time t in (a,d] attains its supremum in the right limit at t. Boundary conventions change no integral. The aggregate support volume in (1) is V(s)=sqrt(b)H(s). Consequently

\[
 Z_{[a,d]}:=\int\sum_J\sup_{a\le s\le d} k_s(y,J)dH(y)
 =1+\int_a^d\frac{dV(s)}{V(s)}
 =\boxed{1+\log\frac{H(d)}{H(a)}}.                 \tag{2}
\]

H is continuous and locally absolutely continuous for b>=2. This includes b=2, where H(s)=s+2(s-c)_++(s-2c)_+. Its derivative at a knot need not exist, but log H is absolutely continuous on [a,d], so there is no additional jump term. Where derivatives exist,

\[
 B_b(c,s)=\partial_s\log H(s)
 =(b-1)E[(s-cK)^{-1}\mid S=s].                    \tag{3}
\]

The distributional kernel derivative has total variation 2B_b, from loss of height on the old support and gain on the new shell. Only the new shell enters (2). Replacing it by an asserted TV bound B_b would be wrong.

For the positive conditional convolution K_c^s on H_b times the label group,

\[
 \|\sup_{s\in[a,d]}|K_c^s F|\|_1\le Z_{[a,d]}\|F\|_1.
 \tag{4}
\]

This constant is sharp for this maximal convolution: take an approximate identity concentrated at the neutral label and the origin. Convolution converges to each k_s(-X,J) almost everywhere for a countable dense set of s including endpoints. Fatou and then monotone convergence over finite subsets give the lower bound Z, while (4) gives the upper bound. Density in s recovers the same envelope by the one-sided shell description. This is an operator-norm statement, not an efficient-halo construction.

In particular the loss cannot be uniformly bounded when a decreases to zero. If d<c, Z=1+(b-1)log(d/a). No full-(0,infinity) strong-L1 envelope follows. Nor is (2) a c-maximal bound.

## Uniform macroscopic rate and critical windows

If c and s/b range over fixed positive compact intervals, then uniformly as b tends to infinity,

\[
 B_b(c,s)=\gamma(c,s/b)+O(b^{-1}),                 \tag{5}
\]

where gamma>0 uniquely solves

\[
 s/b=1/\gamma+c/(1+e^{\gamma c}).
\]

Here is a direct justification with all constants uniform on these compacts. Exponential tilt by theta=1-gamma makes the increment mean s/b; gamma stays in a positive compact interval because 1<=1+t/(1+e^t)<=mu. The tilted variance is bounded above and away from zero, and its characteristic function is bounded by (1+t^2/gamma^2)^(-1/2). The centered logarithm has expansion -v t^2/2+O(|t|^3) on a common neighborhood of zero, with all required higher moments bounded uniformly. Fourier inversion gives the tilted density at its mean (2pi bv)^(-1/2)+O(b^(-1)). In its first derivative the Gaussian term is odd and integrates to zero; the local error is bounded by C b integral |t|^4 exp(-k b t^2)dt=O(b^(-3/2)), and the remaining frequencies are exponentially small. Thus the tilted logarithmic density derivative is O(1/b). Undoing the fixed tilt gives f'/f=gamma-1+O(1/b). Since B=1+f'/f, (5) follows. One differentiates the density with the tilt fixed at the evaluation point, not the implicit gamma(s).

For c=mu, mu=1+exp(-mu), put q_b=b mu-sqrt(mu b)z_b with barPhi(z_b)=sqrt(mu/b)z_b. The untilted variance is mu, z_b=O(sqrt(log b)), and the inverse mean map gives gamma=1+O(sqrt(log b/b)) throughout q_b+O(log b). Hence a window of length ell_b=O(log b) satisfies

\[
 Z=1+\ell_b\{1+O(\sqrt{\log b/b})\}.
 \tag{6}
\]

For ell_b=w log b with fixed w>0, log Z=log log b+log w+o(1). The exact formula holds in every finite dimension; (5)-(6) are proved asymptotics, not extrapolations from numerical samples.

## Global physical coordinates and the scope of the transport

For v_i=-log|x_i| in R, let V=sum v_i and X=v-(V/b)1. On each orthant

\[
 dx=e^{-V}\,dv=\frac1{\sqrt b}e^{-V}dV\,dH(X).
\]

When uniform Haar probability is used for orthant labels, the physical modular of F becomes

\[
 \int\Phi_b(f)dx=\frac{2^b}{\sqrt b}
 \int_{\mathbb R}e^{-V}\int\Phi_b(F(V,X,\sigma))dH(X)d\lambda(\sigma)dV.
 \tag{7}
\]

There is no V^(b-1) factor here. That factor belongs to normalized simplex angles, and must not be inserted a second time in this translation-invariant coordinate system. The conditional physical displacement is exactly X_src=X_rec+P_b Z.

The separate Fubini bridge note applies (4) to every fixed source fibre w before taking the receiver supremum. Its dominating source-cost density is covariant under V-translation induced by physical dilation. It thereby removes the fixed-band radial graph trace for nonnegative, band-integrated source cost. It does not prove an action-independent active Bregman baseline, capture the actual successful average inside this band, or compare arbitrary moving physical rectangles.

For product blocks and a product conditioning band the exact kernel envelope masses multiply, by positivity and the tensor product structure. The same Tonelli proof works with vector u,w. In particular fixed r comparable-size critical blocks cost exp(O(r log log N)). No uniform growing-block action theorem is inferred.
