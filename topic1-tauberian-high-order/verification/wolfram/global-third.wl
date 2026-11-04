(* R8 independent scalar certificates for the expanded-window compensation.
   These verify algebraic subclaims; the set-theoretic proof is in geometry.md. *)
ClearAll["Global`*"];
Print["VERSION: ", $Version];
check[name_, expr_] := Print[name, ": ", ToString[expr, InputForm]];

check["expanded_neighbor_radius_violation",
 Reduce[0 < eta < 1/40 && rr > 0 &&
   Max[1,rr] < 1-eta+rr/2 && !(2 eta < rr < 2),
   {eta,rr}, Reals]];

check["thickened_graph_radius_drift_violation",
 Reduce[eta>0 && nn>=1 && 0<h<eta/(8 nn) && 4 nn h>=eta/2,
   {eta,nn,h}, Reals]];

(* log((1+7 eta)/(1-7 eta)) has derivative 14/(1-49 eta^2),
   and vanishes at eta=0. *)
check["regular_edge_log_derivative_violation",
 Reduce[0<eta<1/40 && 14/(1-49 eta^2)>20, {eta}, Reals]];

check["mixed_edge_different_halfbox_violation",
 Reduce[0<eta<1/40 && rr>=3 && pp>0 &&
   pp >= (1-eta)rr && pp/2 < 1+3 eta rr,
   {eta,rr,pp}, Reals]];

check["mixed_edge_same_halfbox_violation",
 Reduce[0<eta<1/40 && rr>=3 && (1-eta)rr<=1,
   {eta,rr}, Reals]];

check["finite_layer_radius_constant", FullSimplify[3 Exp[1/4]<4]];
check["finite_layer_containment_violation",
 Reduce[kk>=0 && 8 kk+2>10(kk+1), {kk}, Reals]];

(* With x=c e in [0,1], the derivative of (1+2x)^n is bounded
   by 2n 3^(n-1). This gives C2(n)=2n c 3^(n-1). *)
check["boundary_protrusion_base_violation",
 Reduce[0<=xx<=1 && (1+2xx>3 || 1+2xx<=0), {xx}, Reals]];
(* Taking the logarithm of the ratio of positive bases gives ll<=0;
   the required exponent nn-1 is nonnegative. *)
check["boundary_protrusion_power_log_violation",
 Reduce[nn>=2 && ll<=0 && (nn-1)ll>0, {nn,ll}, Reals]];

(* Normalize a r_j^n=1: replacing N_j by the private cluster loses
   at most 1-m_j, which is at most d_j. *)
check["boundary_source_replacement_violation",
 Reduce[0<=mm<=near<=1 && dd>=1-mm && near-mm>dd,
   {mm,near,dd}, Reals]];

(* Normalize a=1. The internal multiplicity lower bound is
   SI >= 2 MI-(SI-sigma)-cross. *)
check["internal_shared_mass_violation",
 Reduce[ss>=0 && mm>=0 && sig>=0 && cross>=0 &&
   ss>=2mm-(ss-sig)-cross &&
   ss-mm<sig/2-cross/2, {ss,mm,sig,cross}, Reals]];

check["fixed_D_boundary_constant_reduction_violation",
 Reduce[cc>0 && vv>0 && cc vv<=1 && pp>=0 && cc vv pp>pp,
   {cc,vv,pp}, Reals]];

check["positive_recovery_after_boundary_and_volume_absorption",
 FullSimplify[cc-cc/4-cc/4]];

(* For 0<vv<1, ll=Log[vv]<0; taking logarithms of
   vv^(nn+1)<=vv^3 leaves (nn-2)ll<=0. *)
check["n_ge_2_final_remainder_log_violation",
 Reduce[nn>=2 && ll<0 && (nn-2)ll>0, {nn,ll}, Reals]];

check["target_denominator_identity",
 FullSimplify[1+aa/(bb-aa)-1/(1-aa/bb),
   Assumptions -> bb>aa>0]];

Print["END"];
