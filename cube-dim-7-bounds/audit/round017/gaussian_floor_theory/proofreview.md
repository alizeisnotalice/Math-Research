# R017 independent mathematical audit: Gaussian-proxy paired-budget barrier

2026-10-27. **The named-proxy bounds and paired-certificate barrier are accepted, with the raw-shape conditioning clarification in section1.** The exact barrier is

    K=7^(14/3)pi^3/840>128.

It applies only to the specified moment-based numerator and specified floor-proxy certificates. It is not a lower bound on C7, an upper bound on the actual cube floor, a lower bound on actualI, or exclusion of all Gaussian/Green methods. This audit ran no mathematical numerical or symbolic program, scalar gate, cache functional, phase, kernel or sample. Parent retains proof responsibility.

## 1. Raw shape versus radial-conditioned shape: local wording correction

In proof.md:65, the exact minimum of the 'raw single-shape kernel' is called L_vertex. If raw means Sigma before conditioning away the Gamma variable, that literal pointwise statement is false. With R~Gamma(32), P independent, B=A(P)^(-1), d=sqrt(detB), C=1/(2pi^(7/2)), the raw Gaussian all-time kernel is

    K_(R,P)(x)=2C sqrt(R) d (x^T Bx)^(-3).

Its exact cube minimum is

    ell_(R,P)=2C sqrt(R) d/[max_sigma sigma^T Bsigma]^3.

The defined L_vertex(P)=CD d/[max_sigma sigma^T Bsigma]^3, where D=2E sqrt(R), is its **R-conditional expectation**, not this raw pointwise minimum. For any measurable cube-constant raw-shape minorant m(R,P), m<=ell_(R,P) pointwise, so

    E m<=E ell=E_P L_vertex.

This proves exactly the expectation assertion needed by the barrier. If 'shape' means the already R-conditioned P-kernel, L_vertex is its exact minimum directly. Both interpretations yield the same accepted expectation conclusion once the conditioning is written explicitly. This review does not edit the frozen source or accept the mislabeled pointwise equality as written.

## 2. Two lower-bound mechanisms and their quantifiers

Let P~Dirichlet(1/2)^64, A=sum_v vv^T/P_v, B=A^-1 and d=sqrt(detB). The accepted original common-time representation is

    G(x)=CD E_P[d(x^T Bx)^(-3)].

It is an extended-positive identity; no conditional Gaussian coordinate symmetry is assumed. Define

    L_vertex=CD d/[max_(sigma in{+/-1}^7) sigma^T Bsigma]^3,
    L_trace =CD d/(trB)^3.

For each SPD B, coordinatewise convexity puts the cube maximum of x^T Bx at a vertex. Hence the R-conditioned P-kernel is at least L_vertex for every cube point. Averaging the vertex quadratics over signs gives trB, so max>=trB and 0<L_vertex<=L_trace. The latter comparison alone does not make L_trace a floor of an unaveraged shape kernel.

The trace floor uses the probability law before Jensen. Signed permutations permute the64unoriented body lines and preserve the Dirichlet law, while B conjugates and d stays unchanged. Average the kernel over the finite group. For x!=0, Jensen for u^-3 gives

    avg_Q (x^T Q BQ^T x)^(-3)
       >=[|x|^2 trB/7]^(-3)>=(trB)^(-3)

on[-1,1]^7. At0 the extended kernel is infinite. Finite group averaging commutes with positive expectation without presupposing finiteness. Alternatively, averaging signs alone gives mean sum_j x_j^2 B_jj<=trB and the same final lower inequality. Thus

    G(x)>=g_trace=E L_trace>=g_vertex=E L_vertex

for every cube point. This is a pointwise statement and therefore also supplies the corresponding a.e. statement, but it gives no explicit certified numerical value here. Trace-Jensen is a law-symmetrized expectation floor, not an unaveraged per-shape floor.

## 3. Deterministic envelope and proxy upper expectations

For every positive P summing to1, weighted Cauchy gives

    u^T Au=sum_v(v dot u)^2/P_v >=[sum_v|v dot u|]^2=psi(u)^2.

The accepted sharp psi floor then yields A>=2048I and B<=I/2048. By determinant AM-GM, trB>=7(detB)^(1/7). Consequently

    L_trace<=CD (detB)^(1/14)/343
            <=CD/[343sqrt2048].

The envelope is deterministic for all shapes, applies also to L_vertex and proves both proxy expectations finite. Existing8192cache means are diagnostic values of bounded functionals, not rigorous expectation certificates in the absence of a separate error argument. No statistic or scalar in that diagnostic is needed in this proof.

Let actualM=E detSigma^(-1/14). The exact radial identity is M=D E_P(detB)^(1/14). Taking expectation in the first preceding inequality gives

    g_vertex<=g_trace<=C M/343.

This is an upper bound on the *values of these named lower surrogates*. It is not an upper bound on actual inf_cube G. Indeed actualG>=g_trace points in the opposite direction. The same expectation upper bound applies to raw cube-constant per-shape minorants by section1. It does not automatically apply to the exact minimum of a richer symmetry-averaged shape kernel or to a construction that preserves x-dependence through expectation; their uniform lower surrogates may be larger than this Jensen proxy.

## 4. Universal barrier for the specified paired chain

Set q=7/6 and

    F=1/[105*2^(25/6)*pi^(13/12)],
    C=1/(2pi^(7/2)), H=16pi^3/15.

The accepted full-boundary argument has I<=F M^q and F=C^q H/896. Suppose a certificate chooses a proved finite M0>=M>0, uses I_cert=F M0^q, and obtains 0<g_cert<=g_trace (or one of the smaller permitted expectation surrogates). Then

    Cost=7^q F M0^q/g_cert^q
        >=7^q F(343/C)^q (M0/M)^q
        >=7^q F(343/C)^q.

Using the boundary form of F removes every fractional pi power and isolates the exact constant:

    7^q F(343/C)^q
      =(7*343)^q H/896
      =7^(14/3)pi^3/840=K.

The factors343=7^3 and896come from determinant/trace AM-GM and the complete cube boundary respectively. There is no dependence on the sharp floor c in K; the sharp floor was only used above to supply the bounded diagnostic envelope. The moment numerator and proxy ceiling must involve the same actualM; substituting unrelated moment surrogates without proving the stated relations is not covered.

The exact analytic comparison to128 is valid without a numerical scalar run. Since7^(2/3)>3 follows by cubing, and pi>3 follows from an inscribed regular hexagon,

    K>7^4*3*27/840=194481/840>128.

The final rational inequality is194481>107520. These are symbolic proof steps in the source, not new numerical-tool output. This barrier persists even at M0=M. It is a lower bound on a reported *upper-certificate expression*. As I_cert is an upper bound, not actualI, a cost lower bound cannot be promoted to a C7lower bound.

## 5. R016 specialization and the precise numerator requirement

The R016 budget J0=4018421141119/10^18 was constructed by rounding upward an expression F M_RB^q with a proved moment bound M_RB>=M. Therefore defining M0=(J0/F)^(1/q) indeed yields M0>=M_RB>=M. This follows from R016's moment provenance; the implication would be invalid for an arbitrary unrelated Iupper merely satisfying I<=J0.

For this same numerator, let U_proxy=C M0/343 and g_req=7(J0/128)^(1/q). Algebra gives

    U_proxy/g_req=(128/K)^(1/q)<1.

Thus no improved estimate of either named proxy expectation can certify a floor large enough for this same numerator chain to reach128. The already reported finite-cache means do not enter this exclusion and need not be recomputed. The four-scalar gate in registration.md is only a finite calibration of a theorem whose strict comparison is already proved analytically. This independent audit did not run that gate or rely on its future outputs.

## 6. Strict accepted range and limits

Accepted: uniform vertex floor after the correctly identified conditioning; trace floor after law symmetrization and Jensen; hierarchy and deterministic envelope; expectation ceiling C M/343; exact K and strict K>128; extension to raw cube-constant shape minorants after conditional expectation; R016 specialization by its actual moment-built numerator. The raw-floor wording in section1 must accompany acceptance.

Not implied: actual cube floor<=U_proxy, actualI=F M0^q, actualI lower bounds, C7>128, failure of the original Green model, failure of all Gaussian/x-dependent floor methods, or a rigorous expectation value inferred from the8192cache diagnostic. A different numerator inequality or a richer expectation-before-infimum construction falls outside the proved chain unless separately shown to obey its hypotheses. Globalg=1/320000 remains a separate issue; R016's conditional operator cost cannot be promoted here.

All source documents are read-only. Source manifests and declared dependency bytes were hash-checked separately as metadata. No mathematical numerical or symbolic software, MATLAB, MC, phase, kernel, scalar gate, new cache functional, commit or push was used. Parent retains integration responsibility; this pure audit is complete and stops.
