# R061 disjoint layers: independent full-chain review

2026-10-28. **ACCEPT.** For all sufficiently large integer n, the final frozen proof establishes

    L=log n, x=47L/50, t=n/x,
    C_n<=F(A_n(t),B_n(t)-(429/2000)n)
       <n log n-n log log n-3n/20.

This is an original Lebesgue centred axis-parallel cube maximal-operator bound for arbitrary L1 input and all positive common radii. It has no explicit starting dimension. No R056 or R063 threshold, n-dependent layer refinement, ideal-model optimum, or lower bound for C_n is imported. The main coordinator retains final integration/proof responsibility. No blocking gap or requested author correction was found.

Final audited proof: `/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round061/disjoint_layers/PROOF.md`, actual SHA-256 `0a1e2bc6eb829714f62b9a0cdb5f52f651fa7d49fd54a084e89b2ff3751af716` (345 lines). Only this independent clone's `audit/round061/disjoint_layers/` was written. No author/main file, Git state, production run, grid, n=7 work or MATLAB was used or changed.

## 1. Read dependencies and division of responsibility

The full final proof was read and its complete operator chain rebuilt. The following were also read: the independent general direct-sum proof/counterexamples, the R062 independent arithmetic audit's code/registration/receipt and certificate interface, and the accepted R056/local-J statements used by the new geometry. Actual locks are recorded in `microcheck.json`:

| Input | SHA-256 |
|---|---|
| R061 full source | `0a1e2bc6eb829714f62b9a0cdb5f52f651fa7d49fd54a084e89b2ff3751af716` |
| General direct-sum lemma | `8b8168370146670889ad5568fc3b20c48df5e020fbe64da8ffac29252dbf6f3f` |
| R062 fixed 128-layer certificate | `b889e35ef909bd979937752c67410bf3913eb76cdf7dc6c720556e093d4c0711` |
| R062 independent arithmetic receipt | `77cc665be3192d3532aabdf9879e827a49748c20e072b79eb0b422ec6582cf54` |
| Expanded local J lemma | `97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec` |
| Accepted R056 full mechanism | `6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a` |

The independent arithmetic code reconstructs all fixed row fractions and their 48-term log/rounding directions, enclosing them separately by 40-term intervals. Its actual code and registration hashes match its receipt. That conditional arithmetic input does not itself prove geometry or an operator bound; those missing links are checked below. This audit did not execute either author reconstruction, the numerical audit, its 512-layer ideal-model portion, any 128/129-row log baseline or a producer helper.

## 2. One deficit measure, not a sum of item weights

For one actual mu<=lambda, set nu=lambda-mu. For each layer, expansion gives

    E_lambda,lambda,j-E_mu,mu,j
      =E_nu_Aj,lambda_Pj+E_mu_Aj,nu_Pj
      <=M_j k_j[nu(A_j)+nu(P_j)].

The actual row bound gives E_mu,mu,j<=tau*mu(A_j)<=tau*lambda(A_j), and the ambient row lower gives E_lambda,lambda,j>=lambda(A_j)r_j k_j. Therefore

    nu(A_j union P_j)
       >=lambda(A_j)(r_j-tau/k_j)/M_j.

Here Aj intersect Pj must be null, as explicitly required in the source; it is not a conclusion merely from disjointness of different Bj. In the application **all 2J projected point supports** are disjoint. Summing their missing masses is therefore valid with their own M_j,k_j. Disjointness only in product space would not suffice, as the separate counterexample correctly demonstrates.

This identity does not assign a whole original item to one layer. If one increment crosses several bands, its restricted measures are on disjoint point sets and their sum is at most its original weight. A target item meeting several Aj likewise has its measure portions integrated separately. Every row inequality uses the same actual mu/history and the same global tau. No independently optimized histories or repeated whole-item payments are added.

## 3. Genuine disjoint geometry and finite-h active reference

The fixed partition step is 3/1600. Source intervals partition [18/25,24/25] by exact consecutive endpoints, with shared endpoints assigned once. Target layer j is trimmed by delta_r=2^-30 at each end of its translated interval. Its radial width is step-2delta_r>0; adjacent target intervals have gap 2delta_r. All target radii lie in (1,31/25), hence below 2. For sufficiently large t, all source radii exceed 4. Thus target-target, source-source and every target-source pair are disjoint in point space. Source endpoint spheres are lambda-null and also mu-null by domination. One common terminal-radius normalization is used for all 128 layers.

The matching interval [r-d-h,r-d+h] is strictly inside that target's own source interval when h<delta_r/4. Its ordinary width is exactly 2h, with log width in [2h/qhi,j,2h/qlo,j]. For every active matched source,

    R=q+d-h>=R_j+delta_r-2h>=R_j=qlo,j+d>=1.

The first R_j is exactly one. At the unmatched low edge R can be 1-h<1; the proof does not apply a local R>=1 lemma there. Its b_j extension remains well defined because R>1/2, a0>0 and ell_R>0 eventually, while the unmatched column is zero. This is the necessary finite-h distinction, not blind replacement by a limiting radius.

On active matches q>1/2, R in [r-2h,r], rho>=4r, a0=(1-zeta*t)q>0 and g=d-h+zeta*rho lies in (0,d] once zeta*t<=h/2. Here d=7/25<1/3 and g<R/2. Every guard of the independently proved expanded local J lemma holds, uniformly across the finite layer system. The correct exact-mean baseline is qb=r-d+h with R=r.

## 4. Ambient target mass, pins and calibrated inverse rows

At each full target face, n-1 free coordinates have the accepted exact EJ=1 baseline, log mean/variance estimates, and one deterministic pinned coordinate. Their log-sum mean/variance are O(L^2/n), with the pin term O(L/n), uniformly over all target radii. The fixed log tolerance e=2^-20 is legitimate. Both endpoint counts retain the exact means n+(n-1)a0(qpm)/r and variance <=n/4; no free law is assigned to the pin.

Here 4d>1, so K_w<=n<=4nd is genuinely automatic, unlike the R059 wide tuple. The original target cap complement is exponentially small **per signed face**. These baseline/two-count/cap probabilities all tend to zero uniformly. Thus for any fixed epsA their union retains at least 1-epsA of each full face, giving lambda(A_j)>=(1-epsA)n log(rhi,j/rlo,j)>=m_j n. There is no factor J probability loss: these are separate deterministic mass estimates on disjoint radial sets, and the maximum of finitely many uniform thresholds suffices.

The a.e. derivative bound is integrated using absolute continuity across capture-onset kinks. Since h=2^-48/(Ld), the total log change is <=1280*2^-48. Adding the baseline tolerance yields the common exact P0. Individual J converges uniformly to one and separately enters [1-e,1+e]; it is not used as an nth-power product replacement.

Endpoint source-sign counts and ell_* are monotone. Lower/upper cross-endpoint calibration retains the 1/n pin correction and exact difference 2h(1-zeta*t)(d-h)/[r(r-2h)]<=2hd. With 3h+1/n<=e/2, normalized T_N/(2n ell_*) lies in [1-e,1+e] over the entire matching interval. Consequently the **full signed source-face sum**

    U=(product J)/(2ell_*) sum_i N_i/J_i

lies between the stated common Ulo*n and Uhi*n. Its pin factor has not been discarded. All these full-face statements concern the ambient envelope, not the unknown actual increment distribution.

## 5. Layerwise inverse retention, cap and tau

For layer j, active R>=R_j and g<=d give Qpin<=phi_j<=phibar_j; also R/q<=R_j/qlo,j for finite h. The exact full-source weighted phi moments, divided by J_i>=1-e, bound each free inverse-capture coordinate. On a **fixed captured full source face**, the n-1 free coordinates are independent and the pin is deterministic. Thus the mean is <=phibar_j+x*mubar_j and variance <=x*V_j.

Because Delta_j=chi_j-mubar_j>0, x>=2phibar_j/Delta_j gives a mean gap >=x Delta_j/2. The second displayed condition x>=4V_j/(epsRet Delta_j^2) then bounds that face's tail by epsRet. Uniform face-probability inequalities may be mixed; no within-face variance is falsely equated to total mixture variance. This is done before D restriction and soft conditioning. The bounded source b_j equals theta_j/G on the retained event and is nonnegative elsewhere, so the complete soft row has the stated retention and width factors.

The local cap width is finite and can explicitly be bounded by

    Csrc,j=(qhi,j-qlo,j)/qlo,j,
    n log(qhi,j/qlo,j)<=Csrc,j n.

The original cap source exponent is (1-1/n)eta_j x, with eta_j=(1-qhi,j)/(2qhi,j). Thus b_j<=1 pays the lost row by its unweighted full mass, and division by k_j produces at most Csrc,j/h times exp(-[(1-1/n)eta_j-chi_j]x). This tends to zero when eta_j>chi_j. There is no unproved replacement by a width-independent old cap estimate.

The one global tau=sqrt[A0(B0-gamma n)] is positive eventually: the original cap bounds give B0~nL, while gamma is fixed below one. Its ratio to each k_j is bounded by (sqrt L/h) exp(-[sigma/2-(1-sigma)chi_j]L). Every such exponent is strictly positive. The finite certificate's minimum exact cap/tau gaps are respectively 6379/3000000 and 5191/25000000; retention gaps are at least 1/10000. Hence cap<=epsS k_j and tau<=epsT k_j eventually in every layer. No R056/R063 finite domain is used.

The refined pointwise target face-count/overlap coupling gives column <=2k_j when nonzero. Such sources are active and R>=R_j>=1; inactive columns vanish. The upper row is 2Uhi*k_j/qlo,j, so M_j=2Uhi/qlo,j>2 is valid. The final Borel kernels are explicitly supported in Aj x Pj, bounded by one, with b_j defined on the complete source band before deleting its cap. Their row excess after cap and threshold payment is b_j*k_j, with b_j>0.

## 6. Actual finite histories and the ca scaling

At each load point there is one normalization, one actual uncaptured measure and one global tau. The full finite envelope dominates the actual measure. Restricting both to D gives the exact ca*lambda envelope on the new supports, ca=alpha^-n>=1. Missing mass of the full envelope is at least the missing mass on this restriction; actual mass outside D is not silently removed from the original load.

An actual target point belongs to an item whose uncaptured-centre inequality forces item radius <3. For an actual source increment point z=y/(mathfrak a+1), original positive load support and increment inclusion give ||z||<r_i, so its item radius is >4. Every included source therefore precedes the same actual target, regardless of how either item intersects radial bands. For mathfrak a>=8n^4, eventually t<n, beta<2 and the centre guard is n^-3-(beta+2)/(mathfrak a+1)>0. Point capture implies true centre capture in the original eligibility test. Integrating only included increment pieces and multiplying by b_j<=1 preserves that target's tau budget.

The energy identity is applied to the same restricted actual measure in every layer. Ambient energy scales by ca^2 and degree/target mass by ca. Its exact local finite-a lower bound is

    lambda(A_j)[ca*r_j-tau/k_j]/M_j
       >=lambda(A_j)[r_j-tau/k_j]/M_j.

There is no residual ca^-1 loss. Applying the direct sum to the truly disjoint projected supports gives total missing cap mass >=n sum gamma_j for every original finite history/dilation after the common threshold. Whole item weights are never assigned afresh to each layer.

Adding the original captured tau and terminal self-weight and retaining the original cover factor gives the stated uniform finite-a bound. At fixed n,t,tau only these numerical inequalities are sent to mathfrak a->infinity. Arbitrary strict-level positive witnesses of the same general L1 input, compact finite covers, their selected dilated witnesses and Lebesgue inner regularity remain unchanged. No selector convergence, independent history gluing or scale/input restriction enters. This proves the original-operator root inequality.

## 7. Common threshold and complete scalar cost

All data and J=128 are fixed independently of n. Geometry improves as h tends to zero and t tends to infinity; zeta*t/h is O(n^-2). Baseline/count failure probabilities and individual-J errors tend to zero uniformly. Each fixed positive Delta_j gives finite inverse-retention thresholds. Each positive cap/tau rate dominates its polynomial prefactor. The finite maximum of these conditions, the B/root positivity condition and the final scalar error condition yields one common integer N. This does not prove a numerical N or uniformity under a varying layer record/J.

The chosen gamma=429/2000 is strictly below the certified total. At sigma=3/50, the complete B price is sigma/(1-sigma)=3/47. Therefore

    gamma-3/47=14163/94000,
    gamma-3/47-3/20=63/94000>0.

The complete positive A/cross/constant error divided by n is bounded by exp(-sigma L)+2sqrt L exp(-sigma L/2)+exp(-L), which tends to zero. Including its strict comparison with 63/94000 in the same eventual condition list proves the displayed strict -3n/20. The result is not inferred from sum gamma alone or from an ideal model.

## 8. New independent compressed budget verification

`microcheck.py` was actually run, PID 93558; 11 new exact Fraction certificates PASS. It does not reconstruct the author's 128 individual logs or import any producer helper. It statically binds the consecutive ledger endpoints/common profile, checks the exact sum and minimum rate gaps, and supplies a different analytic lower certificate for the total budget using only TWO four-term positive-log sums.

Here is the new aggregate proof. Put s=3/1600, delta=2^-30, and split the 128 layers into two groups of 64. After writing

    b_j/M_j=(1-epsRet)(Ulo/Uhi)(qlo,j/qhi,j)-e*qlo,j/Uhi,

use qlo,j/Uhi<1 and the increasing q/(q+s). A lower coefficient for each group's (1-epsA) times logarithmic mass is

    c(q)=(1-e)[(1-epsRet)(Ulo/Uhi)q/(q+s)-e],
    q=18/25 or 21/25.

The first 64 untrimmed radial intervals cover [1,28/25], and the next 64 cover [28/25,31/25]. Each group loses only length 128delta through trimming. Since 1/r<=1, its log mass is at least log(28/25)-128delta or log(31/28)-128delta. Four positive atanh terms at u=3/53 and u=3/59 give strict lower logarithms for these two ratios.

The per-row 48-term log tail is <2^-970, since its atanh argument is <=s/2=3/3200<2^-10. Downward rounding loses <2^-200 more, so each L_j is below the true log by <2^-199. The coefficient multiplying each L_j in gamma_j is positive and <1. The final gamma rounding loses <2^-200. Conservatively the total loss from all 128 rows is <128*2^-198=2^-191.

Combining these bounds gives a two-group rational lower sum strictly greater than

    429/2000+1/40000.

Thus even this coarser, independent aggregate budget exceeds the chosen gamma with strict margin, without reproducing the 128-row baseline. The static ledger's exact total also matches the already independently audited value and exceeds gamma. The numerical audit's per-row logarithm/rounding validity remains a distinct accepted arithmetic input, not a claim made by this compressed check.

A small preliminary Fraction probe developed this two-group bound before the registered run; its optional floating display was diagnostic only. The formal receipt contains exact fractions only. A separate read-only identity comparison verified that the numerical audit's current script and registration match the hashes stored in its actual receipt; that program was not executed here.

## 9. Final status and limits

All required projected supports are genuinely disjoint; cross-band increments are paid pointwise against one missing measure; active-R, per-face inverse probability, cap width/rates, common tau, finite-a scaling and original witnesses all close. There are no outstanding mathematical corrections. The final lock is `0a1e2bc6eb829714f62b9a0cdb5f52f651fa7d49fd54a084e89b2ff3751af716`; its earlier two TeX thin-space substitutions do not alter mathematics. Only this fixed finite construction is accepted. Later ideal frontiers, adaptive refinements and explicit-start proposals remain outside this review.
