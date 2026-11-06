# SUPERSEDED DRAFT — r028 formal report proves the stronger arbitrary-input two-ray theorem

Authoritative: /Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal/strong-maximal-endpoint/rounds/r028/report.md. The old coordinate-even argument below is only a subcase.

# r028 fixed global-ratio anchored endpoint — proof draft

For fixed r in[0,1], let k_r(v)=(1+1_{0<v<r})/(1+r),0<v<1.
A_i,r f(x,z)=∫0^1k_r(v)f(vx,z)dv onpositiveorthant; T_r is sameformulaononeaxis.
Fix boundedcompactinput. Fori unused coordinate, everyintermediatefunctionbounded, compactinunusedcoordinates, hasfinitepositivelevelsets andisdominatedbyHproductmin(1,2R/x_used). Thusallfibresincoordinatei areL1andbounded; rearraystepbelowlegal.

Lettildef beeachfiber'sdecreasingrearrangementin x_i, withotherz fixed. Becauseweight x^-1 k_r(t/x)1_{0<t<x} decreasingt, A_i,r f(x,z)≤A_i,r tildef(x,z) pointwise. tildef globallyequimeasurablewithf. Thelatteroutputdecreasesinx. Atthresholds>0 itshaloD(z) lengthd_z isfinitea.e andtotalDfinitebyabovebound(sameupperHproductindicatorcurrentcoordinate). For0<eps<1 choosea_z=(1−eps)d_z (zeroifd_z0), andD_eps=∫a_z. AtsuchinteriorpointA_i,r tildef(a_z,z)>s wherea_z>0 (endpointpossiblyrootnotused). Multiplybyslicelength:
s a_z <∫tildef(t,z) k_r(t/a_z)1t<a_z dt.
Integratez ->s D_eps ≤∫tildef W, W=(1_{t<a_z}+1_{t<r a_z})/(1+r).
Wglobalhaslowlevel1/(1+r),totalpositivesupportD_eps (r0or1degeneratehandled), highlevel2/(1+r)onmeasure rD_eps. ThusWglobaldecreasingrearrangement equalsk_r(y/D_eps)1y<D_eps. HardyLittlewoodglobal⇒
sD_eps≤1/(1+r)[∫0^D_eps f*(y)dy+∫0^rD_eps f*(y)dy]=D_eps T_r f*(D_eps).
Leteps↓0. T_r f*continuousatpositiveDsinceboundedinput andfiniteintegrals, hence s≤T_r f*(D). SinceT_r f*decreasing, itsweaklevel{≥s}contains[0,D]. To getstrictDcomparisonavoidplateau: repeatargumentatthresholds'>s, then T_r f*(D(s'))≥s'>s, soD(s')≤|{T_r f*>s}|. Let s'↓s, D(s')↑D(s) bystrictlevelunion. ThusD(A_i,r f;s)≤D(A_i,r tildef;s)≤D(T_r f*;s). Therefore(A_i,r f)*≤T_r f* a.e.

Iteratepositivetransforms, globallyrearrangingeachphysicalaxisandusingorderpreservingT_r, toobtain(A_n,rn...A_1,r1 f)*≤T_rn...T_r1 f*. Linearproductphysicaloperator isactualfixedratioanchoraverageofcoordinatewise-even extensionf_even toR^n: interval[−r_i x_i,x_i] forx_i>0 andreflectionfornegative; denominatorssame. Everyevenf correspondspositiveorthantf, distributionandmodularbothmultiply2^n, cancelsratio.

For fixedvectorr andcapA>1 defineQ(n,r,A)=sup{1/J(h):h:[0,∞)→[0,A] nondecreasing, E h(S_r)>1},J=∫e^-tΦ_n(h)dt.
BecauseTproducts g(x)=E g(xe^-S) andg=f*decreases, anyxstrictlyinsideitsunit-halo satisfiesE f*(xe^-S)>1. h(t)=f*(xe^-t) isboundedA,nondecreasing, J(h)=x^-1∫0^xΦ(f*)≤modular(f)/x. Henceoutputweakratio≤Q(n,r,A). Physicalinputtruncationmin(f,M)1box andmonotoneconvergenceextendgeneralPhi-integrableevenf(capAalreadybounded,truncateboxonly).

Lower: takeanyfeasibleh boundedA, extend0below0. Actualevenf_L(x)=h(Σ−logabsx_i−L)1[-1,1]^n. ForxinsidecubeandSx≥L, chosenfixedranchoraverage=E h(Sx−L+S_r)≥E h(S_r)>1. Gammaexactvolumeandcostdeepasymptgiveweakratio≥1/J(h). Thereforebestactualfixedratioanchoredconstant EXACTQ(n,r,A). SupvectorsgivescapQ(n,A); allowA∞ viafinitecaps/truncatehsoequalsQ_n. NoM_Supper,noinput-dependentglobalratio,nooutputpointwisechoice,noarbitrarynoneveninputextensionclaim.

Criticalstrictplateaufixabove(s'>s)mustappearformalreport. ScopeincludesALLmicroblock+multiple-receiverprofilesatfixedcommonc, notonlysinglecone. Thusstopallfixedglobalrgeometrysearchifaccepted.
