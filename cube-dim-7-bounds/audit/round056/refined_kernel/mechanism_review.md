# R056 independent refined-kernel mechanism review

Conclusion: **ACCEPT the assigned local mechanisms**, at the frozen source below and subject to the separately audited local-J inputs and finite guards. This is not a separate final acceptance of the explicit operator threshold or final Cn inequality; the full-chain/finite-threshold auditor owns that conclusion. No new Cn claim is activated by this scoped report alone.

Frozen author source:

    /Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md
    SHA256 6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a

The author supplied the initial 5d639f31 source, followed by this final frozen revision. The final changes were inspected: its target-cap probability is stated per face without an invalid averaged factor 1/2, and it explicitly defines the full source band and b before D restriction, the final kernel, degrees, disjoint sets and actual source-radius relation. The assigned local mechanisms and constants are unchanged. No author file was modified. This audit writes only the independent clone's `audit/round056/refined_kernel/` directory and does not rerun the local-J task, the full 29-gate checker, an old grid, a production scan, MATLAB, or Git operations.

## 1. Required inputs, distinct probability spaces and scoped checks

The mechanism uses d=3/10, q in [31/40,37/40], target r in the trimmed interval given by the source, h=2^-24/(Ld), x=93L/100, and theta=exp(-7x/200). The local coordinate task supplies its stated J/product bounds at d<=1/3 under R>=1, g<=R/2 and rho>=4R; this audit does not accept a domain extension solely by citing an old d<=1/16 result.

There are three different uses of probability, which the frozen source keeps separate:

1. full target-face coordinate laws establish baseline log-product and two endpoint count events;
2. for a **fixed good target**, each captured full source face has a product inverse-capture law;
3. the original D cap is imposed only after the inverse-face retention estimate, through a deterministic unweighted-mass loss.

No independence of coordinates after D conditioning is assumed. No equality between a mixture's variance and its within-face variances is used.

## 2. Exact target column factor 2

Fix a source z. Because |z_i|<=rho<kappa*rho and kappa*rho>r, put a_i=kappa*rho-|z_i|>0. A target coordinate capture interval contains zero. Its normalized side fraction is

    ell_i(r)=min(1,(1+a_i/r)/2).

If a_i>=r, ell_i=1 and at most two target faces are captured. If 0<a_i<r, exactly one target face is captured, and ell_i>1/2. At equality a_i=r the strict-boundary face convention can only reduce the count. Hence pointwise

    N_target,i/(2ell_i)<=1.

The exact target shell coefficient is

    (product ell_i)*sum_i N_target,i/(2ell_i)
      <=n*product ell_i.

This is the correct factor n; the old independent deleted-coordinate bound 2n is unnecessarily weak here. Monotonicity of each ell_i as r increases from R and matched log-width at most 2h/R give

    C_raw(z)<=2nh*G(z)/R.

With b=min(1,theta/G), bG<=theta and R>1, the final column is at most 2k. This estimate holds before any source concentration and after every target-set restriction by monotonicity. It does not use an average face count or a probabilistic replacement for the column maximum.

## 3. Source pin-count calibration and the full row coefficient

Distinguish the source pinned-face count in the row from the target-face count in the preceding column calculation. Here

    a0(q)=(1-zeta*t)q,
    N_i(q)=1+1_{|w_i|<a0(q)},
    C(q)=sum_i N_i(q),
    ell_pin(q)=[R(q)+a0(q)]/[2R(q)].

On a fixed maximal target face, its pinned coordinate has N_i=1 since a0(q)<r. The n-1 remaining indicators are independent Bernoulli with probability a0(q)/r, so

    E C(q)=n+(n-1)a0(q)/r,
    Var C(q)<=n/4.

Chebyshev at tolerance e*n/4 gives probability at most 4/(e^2*n) for each endpoint. The n-1 correction is retained in the means and later contributes the 1/n drift.

Both C(q) and ell_pin(q) increase with q. At q+=r-d+h, R=r; at q-=q+-2h, R=r-2h. On the two endpoint count events, monotone cross-endpoint bounds control every intermediate q. The lower cross-mean ratio loses at most 2h+1/n. For the upper cross-mean ratio the exact difference is

    a_+/r-a_-/(r-2h)
       =2h(1-zeta*t)(d-h)/[r(r-2h)]<=2hd.

The random count error divided by 2n*ell_pin is at most e/4 because ell_pin>=1/2. Therefore the stated gate 3h+1/n<=e/2 yields the valid conservative uniform enclosure

    1-e<=C(q)/(2n*ell_pin(q))<=1+e.

This does not sum independent worst positions. The same fixed target w supplies both endpoint events and the entire matching interval.

The complete inverse-source shell coefficient is exactly

    U=(product J)/(2ell_pin)*sum_i N_i/J_i.

Combining the fixed product enclosure P0<=product J<=1/P0, the individual J enclosure 1-e<=J_i<=1+e, and the normalized count enclosure gives precisely

    Ulo=P0(1-e)/(1+e),
    Uhi=(1+e)/[P0(1-e)],
    Ulo*n<=U<=Uhi*n.

The pin factor is not dropped: it is canceled by the concentration of its corresponding full face count. There is no additional independence assumption between the count and product, because all inequalities are pointwise on the same good-target set.

## 4. Per-source-face softclip Chebyshev, then a legitimate mixture

For a fixed captured source face, z_j=+-rho is pinned. Its contribution to Phi=-logG is Qpin=-logell_pin and has variance zero. Its other source coordinates are independent under densities proportional to capture/ell, with normalizers J_i>=1-e. Fixing w in the good target set does not break that product factorization of the full source face.

The free-coordinate full-source integrals are exactly

    E_full[(-logell)/ell]=R*Qpin^2/rho,
    E_full[(-logell)^2/ell]=2R*Qpin^3/(3rho).

They follow by substituting y=1-u/(2R) on the boundary strip u in [0,g]. Dropping the capture indicator only increases these nonnegative numerators. Dividing by the actual inverse-capture normalizer J_i then gives the source's per-face bounds

    E_face Phi<=Qpin+x*R*Qpin^2/[q(1-e)],
    Var_face Phi<=x*2R*Qpin^3/[3q(1-e)].

Independence is used only within this one full source face. Bounding each variance by its second moment is legitimate; the pinned term contributes no variance.

For h<=1/10000, R>=10749/10000 and g<=3/10 give

    Qpin<=log(3583/3083).

Also R/q=1+(d-h)/q<=43/31. With u=500/3083, the third-order Taylor upper bound Q=u-u^2/2+u^3/3 is rigorous for log(1+u). Independent exact checks confirm Q<1/6 and the two stated rate bounds

    [43/(31(1-e))]*Q^2<4/125,
    (2/3)*[43/(31(1-e))]*Q^3<1/300.

At the stated pin-budget guard, every face mean is at least x/400 below the common threshold chi*x. Hence each face tail probability is at most

    (x/300)/(x/400)^2=1600/(3x)<=e.

This bound is uniform over the captured faces and their inverse-row mixture weights. Averaging the **probabilities** therefore preserves e; averaging only the within-face variances would not be enough. The frozen source explicitly uses the former argument.

On G>=theta, soft mass equals theta times inverse mass, and elsewhere its ratio to inverse mass is between zero and theta. Thus the full soft shell coefficient is between `(1-e)Ulo*n*theta` and `Uhi*n*theta`. This interpretation accounts for the soft clipped tail without equating two different probability laws.

## 5. Full matching width and deterministic post-cap payment

The trimmed target band and h<=2^-20 ensure q-=r-d-h>=qmin and q+=r-d+h<=qmax. The ordinary matched width is therefore exactly 2h, not the old one-sided h lower bound. Its logarithmic width lies in [2h/qmax,2h/qmin]. The full soft row bounds follow with factors 2Ulo/qmax and 2Uhi/qmin.

The original source cap is then removed, without reusing any conditioned-coordinate independence. Its negative pinned source face is excluded. On a positive face its free fraction is

    ell_t(rho)=1-(1-q)/(2rho).

Since `(1-q)/(2q)>=3/74` and n>=75,

    (n-1)*(1-q)/(2rho)
      >=(74/75)*x*(3/74)=x/25.

The full source-band mass is n log(37/31)<nd. Thus the cap mass is at most nd*exp(-x/25). Because b<=1, this full unweighted mass bounds the lost soft row uniformly for every target. After division by k=nh*exp(-chi*x), its ratio is at most

    (d/h)*exp[-(1/25-chi)x]=(d/h)*exp(-x/200).

This is the valid new cap exponent. The older exp(-Ld/32) would fail to pay the present theta and must not be substituted. The frozen source uses the corrected x/25 estimate.

The independent rate check verifies the exact gaps 1/25-chi=1/200 and chi-4/125-1/2000=1/400. The full-chain auditor owns the Taylor-eight proof that all cap/tau/probability gates hold uniformly for L>=2^20; this report does not rerun that 29-gate computation.

## 6. Degree/actual-history interface

After cap deletion, the row lower budget subtracts epsS*k. Paying tau subtracts epsT*k. These are separate terms: the first is the post-cap loss, the second is the true selection threshold. The upper row remains at most 2Uhi*k/qmin, and exact arithmetic confirms this is below (21/8)k. Columns are at most 2k, so the proposed common degree denominator (21/8)k is valid at this mechanism level.

The kernel still has source weight 0<=b<=1. Therefore its actual original-mu row, on one actual target item, is bounded by that same item's original eligibility tau. The radial bands and matching do not add sources from another history. All target items are below radius 3 and all source items above radius 4 once the stated scale guards hold; the established finite-dilation center guard then gives true capture.

The ambient energy/degrees scale by alpha^-2n and alpha^-n on the disjoint target/source sets, while the actual row is still paid by tau. Applying the same nonnegative-kernel missing-measure identity is therefore legitimate. No source-face concentration statement is applied directly to unknown mu. General-L1 and all-positive-radius validity continue through the original finite witness cover, not through a model-to-operator extrapolation.

No gap was found in this assigned mechanism interface. Acceptance of the final explicit Cn theorem additionally requires the separately audited local-J lemma, all finite guards, the exact scalar remainder and the complete original-operator passage. Those scopes are intentionally not claimed by the word ACCEPT at the start of this local report.
