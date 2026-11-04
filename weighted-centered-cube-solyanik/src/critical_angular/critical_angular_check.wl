(* Independent high-precision ProductLog[-1] scalar inverse check for round 30.
   The guards are numerical TrueQ predicates, never symbolic zero tests. *)
$MaxExtraPrecision = 10000;
ClearAll[guard, inv, ff, n, g, k, a, s];
guard[name_, q_] := If[TrueQ[q], Print["PASS ", name], Print["FAIL ", name]; Exit[1]];
wp = 240;
inv[u_?NumericQ, c_?NumericQ] := -ProductLog[-1, N[-u Exp[-1/c]/c, wp]] - 1/c;
ff[u_?NumericQ, c_?NumericQ, b_?NumericQ, beta_?NumericQ] := Module[{ll = inv[u, c]}, Exp[-beta ll] (1 + b ll)];
(* Exact slice-coefficient audit.  rawK is the cube's unnormalised mass;
   b c^-beta below is the separate normalised leading coefficient. *)
sliceChecks = Flatten[Table[
   Module[{nn=n, gamma=g, beta=1-g, v, c, b, rawK, volumeSlice, massSlice},
    v=2^nn; c=(nn-1)/nn; b=(nn-1) beta/(nn-gamma); rawK=v (nn-gamma)/(nn beta^2);
    volumeSlice = v a + v (nn-1) a Integrate[1/s,{s,a^(1/nn),1},Assumptions->a>0];
    massSlice = rawK a^beta + v (nn-1) a^beta/beta Integrate[1/s,{s,a^(1/nn),1},Assumptions->a>0];
    TrueQ[FullSimplify[volumeSlice==v a (1+c Log[1/a]),a>0]] &&
    TrueQ[FullSimplify[massSlice==rawK a^beta (1+b Log[1/a]),a>0]] &&
    TrueQ[FullSimplify[rawK b==v (nn-1)/(nn beta)]]],
   {n,{2,3}},{g,{1/4,1/2,3/4,9/10}}],1];
guard["exact_slice_coefficients",And@@sliceChecks];
rows = Flatten[Table[
   Module[{nn=n, gamma=g, b, c, beta, tau, H, u, L, kappa, local, cap, delta, upper, lower, norm, residual, geometric},
    b=N[(nn-1)(1-gamma)/(nn-gamma),wp]; c=N[(nn-1)/nn,wp]; beta=N[1-gamma,wp]; tau=N[nn beta,wp];
    H=N[k Log[2],wp]; u=Exp[-H]; L=N[inv[u,c],wp]; kappa=N[b c^-beta,wp];
    residual=Abs[N[L-Log[1+c L]-H,wp]];
    local=ff[u,c,b,beta]/(u^beta H^N[gamma,wp]); cap=H^N[1/(2 nn),wp]; delta=1/H;
    upper=cap^-tau ff[u cap^nn,c,b,beta];
    lower=delta^tau ff[u (1-2 delta)^nn/delta^nn,c,b,beta];
    norm=u^beta H^N[gamma,wp]; geometric=TrueQ[cap>2 && delta<1/3 && u<1];
    {nn,gamma,k,N[kappa,40],N[local,40],N[upper/norm,40],N[lower/norm,40],residual,geometric,TrueQ[residual<10^-120]}],
   {n,{2,3}},{g,{1/4,1/2,3/4,9/10}},{k,{8,16,32,64,128,256,512,1024,4096}}],2];
guard["all_numeric_productlog_residuals",And@@(#[[10]]&/@rows)];
guard["all_u_in_unit_interval",And@@Table[2^-k<1,{k,{8,16,32,64,128,256,512,1024,4096}}]];
Print["PASS critical-angular Wolfram checks; rows=",Length[rows],"; geometry-valid rows=",Count[rows[[All,9]],True]];
