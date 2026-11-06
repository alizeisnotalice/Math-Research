(* r058: frozen-source, radial-compensated receiver response diagnostic.
   Floating finite-window quadrature only; source profile, roots, tau and
   calibrated lambda remain fixed while (d,v,a) vary. *)
root58 = DirectoryName[$InputFileName]; args58 = $ScriptCommandLine;
out58 = If[Length[args58] >= 2, args58[[2]], FileNameJoin[{$TemporaryDirectory, "compute58.json"}]];
env58=Quiet@Check[Environment["STRONG_MAXIMAL_REPO"],$Failed];repo58=If[StringQ[env58]&&env58!="",env58,DirectoryName[root58]];If[!FileExistsQ[FileNameJoin[{repo58,"code","primal52.wl"}]],Print["Set STRONG_MAXIMAL_REPO to project root"];Exit[2]];
Block[{$ScriptCommandLine = {}, $primal52LibraryMode = True},
  Get[FileNameJoin[{repo58, "code", "primal52.wl"}]]
];
seed58 = Import[FileNameJoin[{repo58, "rounds", "r052", "primal52.json"}], "RawJSON"];
seedRow57[nn_] := SelectFirst[seed58["rows"], #["n"] == nn &];

(* Analytic t-derivative of the receiver/source c-score s=partial_c log R.
   This is kept separate from primal52's five-slot jet, whose third slot is
   the density second ratio and whose temporary sums are not s_t. *)
sTBuild58[bb_Integer, cc_Real] := sTBuild58[bb, cc] = With[
  {m = bb - 1, nb = bb + 1, bc = bb, cv = cc, pp = 1./(1. + Exp[cc]),
   lb = Developer`ToPackedArray[Log[N[Binomial[bb, #] & /@ Range[0, bb]]]]},
  Compile[{{tt, _Real}}, Module[{u = 0., uT = 0., uC = 0., uCT = 0., y = 0., ly = 0., z = 0., kv = 0., i = 0, kk = 0., ss = 0.},
    If[tt <= 0., Return[0.]];
    For[i = 1, i <= nb, i++, y = (tt - cv (i - 1))/tt;
      If[y > 0., ly = Log[y]; kv = i - 1.; z = lb[[i]] + m ly; If[z > -700., u += Exp[z]];
        z = lb[[i]] + (m - 1) ly; If[z > -700., uT += Exp[z]; uC += kv Exp[z]];
        z = lb[[i]] + (m - 2) ly; If[z > -700., uCT += kv Exp[z]]]];
    If[u <= 0., Return[0.]]; kk = m uT/(tt u); ss = bc pp - m uC/(tt u);
    (bc pp - ss) kk - m (m - 1) uCT/(tt^2 u)
  ], RuntimeOptions -> "Speed"]
];

(* Solve the common baseline only once, then freeze this source H_0. *)
makeFrozen57[nn_, no_Integer, ord_Integer] := Module[
  {s = seedRow57[nn], cc, hi, ctx, sol, ls, bb, js, tx, tw, ir, tau, cap, asc,
   dat, bd, p, la, roots},
  cc = s["tuned_c"]; hi = s["fixed_window"][[2]];
  ctx = context52[nn, cc, cc, no, ord, hi];
  sol = solve52[nn, cc, cc, no, ord, hi];
  ls = sol["raw_integrals", "log_scale"]; bb = nn/2; js = ctx["j1"];
  tx = ctx["x"]; tw = ctx["w"]; ir = ctx["inner_rule"]; tau = ctx["tau"];
  cap = ctx["cap"]; asc = ctx["aCap"];
  dat = Table[
    p = js[tx[[i]]]; la = ls - tau + p[[5]];
    bd[lev_] := Module[{lo = 0., up = hi, md},
      If[la + js[hi][[5]] < Log[lev], Return[hi]];
      Do[md = (lo + up)/2.; If[la + js[md][[5]] < Log[lev], lo = md, up = md], {48}];
      (lo + up)/2.
    ];
    roots = bd /@ {1., N[nn], asc};
    <|"x" -> tx[[i]], "w" -> tw[[i]], "logAmp" -> la, "roots" -> roots|>,
    {i, Length[tx]}];
  <|"n" -> nn, "b" -> bb, "c" -> cc, "hi" -> hi, "ctx" -> ctx,
    "log_scale" -> ls, "tau" -> tau, "cap" -> cap, "aCap" -> asc,
    "source_j" -> js, "data" -> dat, "inner_rule" -> ir, "baseline_solution" -> sol|>
];

(* G(d,v,a)=E_{c+a,c+a} H0(T1+d+v,T2+d-v).
   In source coordinates x=T+(d+v,d-v), so receiver arguments are
   x1-d-v and x2-d+v.  Fa/Faa are analytic c-jets of receiver law. *)
evalFrozen57[mo_, dd_, vv_, aa_] := Module[
  {nn = mo["n"], bb = mo["b"], cc = mo["c"], hi = mo["hi"], js = mo["source_j"],
   ir = mo["inner_rule"], cap = mo["cap"], jr, ans = {0., 0., 0.}, under = 0,
   row, t1, w1, la, rt, x1, p1, f1, seg, l, r, hh, cj, base, xs, ws,
   t2, p2, f2, amp, gp, sec, j},
  jr = jetBuild52[bb, N[cc + aa]];
  Do[
    t1 = row["x"]; w1 = row["w"]; la = row["logAmp"]; rt = row["roots"];
    x1 = t1 - dd - vv; p1 = jr[x1];
    If[p1[[5]] - x1 < -700., under++]; f1 = expSafe52[p1[[5]] - x1];
    Do[
      {l, r, hh} = seg;
      If[r > l,
        cj = cdfJets52[bb, cc + aa, l - dd + vv, r - dd + vv]; base = w1 hh f1;
        ans[[1]] += base cj[[1]];
        ans[[2]] += base (p1[[2]] cj[[1]] + cj[[2]]);
        ans[[3]] += base (p1[[3]] cj[[1]] + 2 p1[[2]] cj[[2]] + cj[[3]])
      ], {seg, {{rt[[1]], rt[[2]], 1.}, {rt[[3]], hi, cap}}}];
    l = rt[[2]]; r = rt[[3]];
    If[r > l,
      xs = (l + r)/2. + (r - l) ir[[1]]/2.; ws = (r - l) ir[[2]]/2.;
      Do[
        t2 = xs[[j]]; p2 = jr[t2 - dd + vv];
        If[p2[[5]] - (t2 - dd + vv) < -700., under++];
        f2 = expSafe52[p2[[5]] - (t2 - dd + vv)];
        amp = expSafe52[la + js[t2][[5]]]; hh = profile52[amp, nn, cap];
        base = w1 ws[[j]] hh f1 f2; gp = p1[[2]] + p2[[2]];
        (* jet slot 3 is f_cc/f, not a log-density second derivative. *)
        sec = p1[[3]] + p2[[3]] + 2 p1[[2]] p2[[2]];
        ans[[1]] += base; ans[[2]] += base gp; ans[[3]] += base sec,
        {j, Length[xs]}]
    ], {row, mo["data"]}];
  <|"G" -> ans[[1]], "Ga" -> ans[[2]], "Gaa" -> ans[[3]],
    "receiver_log_underflow_samples" -> under|>
];

(* Complete r050 response form evaluated with the frozen source.  The lower
   layer and the soft inverse-curvature part are both retained.  At the
   symmetric base this independently supplies Gd=B(1,k+),
   Gad=B(k+,s+), and Gvv=-B(k-,k-). *)
responseForm58[mo_] := Module[
  {nn = mo["n"], hi = mo["hi"], js = mo["source_j"], ir = mo["inner_rule"],
   cap = mo["cap"], lo = ConstantArray[0., 3], so = ConstantArray[0., 3],
   row, t1, w1, la, rt, p1, t2, p2, fp, kp, km, sp, l, r, xs, ws, amp,
   hh, yy, fac, j},
  Do[
    t1 = row["x"]; w1 = row["w"]; la = row["logAmp"]; rt = row["roots"];
    p1 = js[t1];
    If[rt[[1]] < hi - 10.^-10,
      t2 = rt[[1]]; p2 = js[t2];
      If[p2[[4]] > 0.,
        fp = expSafe52[p1[[5]] + p2[[5]] - t1 - t2]/p2[[4]];
        kp = p1[[4]] + p2[[4]]; km = p1[[4]] - p2[[4]]; sp = p1[[2]] + p2[[2]];
        lo += w1 fp {kp, kp sp, km^2}
      ]
    ];
    l = rt[[2]]; r = rt[[3]];
    If[r > l,
      xs = (l + r)/2. + (r - l) ir[[1]]/2.; ws = (r - l) ir[[2]]/2.;
      Do[
        t2 = xs[[j]]; p2 = js[t2];
        amp = expSafe52[la + p2[[5]]]; hh = profile52[amp, nn, cap]; yy = Log[hh];
        fp = expSafe52[p1[[5]] + p2[[5]] - t1 - t2];
        fac = amp hh/((nn - 1) (nn - 1 + yy) (1 + yy)^(nn - 3));
        kp = p1[[4]] + p2[[4]]; km = p1[[4]] - p2[[4]]; sp = p1[[2]] + p2[[2]];
        so += w1 ws[[j]] fp fac {kp, kp sp, km^2},
        {j, Length[xs]}]
    ],
    {row, mo["data"]}];
  <|"B_1_kplus" -> lo[[1]] + so[[1]], "B_kplus_splus" -> lo[[2]] + so[[2]],
    "B_kminus_kminus" -> lo[[3]] + so[[3]],
    "lower_components" -> <|"B_1_kplus" -> lo[[1]], "B_kplus_splus" -> lo[[2]], "B_kminus_kminus" -> lo[[3]]|>,
    "soft_components" -> <|"B_1_kplus" -> so[[1]], "B_kplus_splus" -> so[[2]], "B_kminus_kminus" -> so[[3]]|>|>
];

(* Direct bilinear-response evaluation of the K numerator.  It avoids forming
   Gavv as a small difference of receiver finite differences. *)
directK58[mo_] := Module[
  {nn = mo["n"], hi = mo["hi"], js = mo["source_j"], ir = mo["inner_rule"], cap = mo["cap"],
   st = sTBuild58[mo["b"], mo["c"]], lo = ConstantArray[0., 5], so = ConstantArray[0., 5],
   row, t1, w1, la, rt, p1, t2, p2, fp, kp, km, sp, stm, l, r, xs, ws, amp, hh, yy, fac, j, all, beta, gavv},
  Do[
    t1 = row["x"]; w1 = row["w"]; la = row["logAmp"]; rt = row["roots"]; p1 = js[t1];
    If[rt[[1]] < hi - 10.^-10,
      t2 = rt[[1]]; p2 = js[t2];
      If[p2[[4]] > 0., fp = expSafe52[p1[[5]] + p2[[5]] - t1 - t2]/p2[[4]];
        kp = p1[[4]] + p2[[4]]; km = p1[[4]] - p2[[4]]; sp = p1[[2]] + p2[[2]]; stm = st[t1] - st[t2];
        lo += w1 fp {kp, kp sp, km^2, km stm, km sp km}]
    ];
    l = rt[[2]]; r = rt[[3]];
    If[r > l, xs = (l + r)/2. + (r - l) ir[[1]]/2.; ws = (r - l) ir[[2]]/2.;
      Do[t2 = xs[[j]]; p2 = js[t2]; amp = expSafe52[la + p2[[5]]]; hh = profile52[amp, nn, cap]; yy = Log[hh];
        fp = expSafe52[p1[[5]] + p2[[5]] - t1 - t2]; fac = amp hh/((nn - 1) (nn - 1 + yy) (1 + yy)^(nn - 3));
        kp = p1[[4]] + p2[[4]]; km = p1[[4]] - p2[[4]]; sp = p1[[2]] + p2[[2]]; stm = st[t1] - st[t2];
        so += w1 ws[[j]] fp fac {kp, kp sp, km^2, km stm, km sp km}, {j, Length[xs]}]
    ], {row, mo["data"]}];
  all = lo + so; beta = all[[2]]/all[[1]]; gavv = -(all[[4]] + all[[5]]);
  <|"Gd_B_1_kplus" -> all[[1]], "Gad_B_kplus_splus" -> all[[2]], "Gvv_minus_B_kminus_kminus" -> -all[[3]],
    "Gavv_direct_minus_B_kminus_stminus_plus_splus_kminus" -> gavv,
    "beta_Gad_over_Gd" -> beta,
    "K_direct_response" -> gavv + beta all[[3]],
    "K_direct_single_form" -> -(all[[4]] + all[[5]] - beta all[[3]]),
    "lower_components" -> AssociationThread[{"B_1_kplus", "B_kplus_splus", "B_kminus_kminus", "B_kminus_stminus", "B_kminus_splus_kminus"}, lo],
    "soft_components" -> AssociationThread[{"B_1_kplus", "B_kplus_splus", "B_kminus_kminus", "B_kminus_stminus", "B_kminus_splus_kminus"}, so]|>
];

sTCheck58[mo_] := Module[{j = mo["source_j"], st = sTBuild58[mo["b"], mo["c"]], t, h = .001, an, fd},
  t = mo["b"] (1. + mo["c"]/(1. + Exp[mo["c"]])); an = st[t]; fd = (j[t + h][[2]] - j[t - h][[2]])/(2 h);
  <|"t_away_from_knots" -> t, "step" -> h, "analytic_s_t" -> an, "centered_fd_s_t" -> fd, "absolute_difference" -> Abs[an - fd]|>
];

response57[nn_, no_Integer, ord_Integer] := Module[
  {mo, z, hp = {.01, .005}, vs = {.5, .25, .125}, vm, dm, gdrows, gd, gad,
   rows, krows, primal, h, gdWard, tplusErr, tplusPass, ga0Small, br, direct},
  mo = makeFrozen57[nn, no, ord]; z = evalFrozen57[mo, 0., 0., 0.]; br = responseForm58[mo];
  direct = <|"Gd_response_B_1_kplus" -> br["B_1_kplus"],
    "Gad_response_B_kplus_splus" -> br["B_kplus_splus"],
    "Gvv_response_minus_B_kminus_kminus" -> -br["B_kminus_kminus"],
    "response_components" -> KeyTake[br, {"lower_components", "soft_components"}]|>;
  rows = Table[
    vm = {evalFrozen57[mo, 0., v, 0.], evalFrozen57[mo, 0., -v, 0.]};
    <|"v" -> v,
      "Gvv" -> (vm[[1, "G"]] + vm[[2, "G"]] - 2 z["G"])/v^2,
      "Gavv" -> (vm[[1, "Ga"]] + vm[[2, "Ga"]] - 2 z["Ga"])/v^2,
      "Gplus" -> vm[[1, "G"]], "Gminus" -> vm[[2, "G"]],
      "Gaplus" -> vm[[1, "Ga"]], "Gaminus" -> vm[[2, "Ga"]]|>, {v, vs}];
  gdrows = Table[
    dm = {evalFrozen57[mo, h, 0., 0.], evalFrozen57[mo, -h, 0., 0.]};
    gd = (dm[[1, "G"]] - dm[[2, "G"]])/(2 h);
    gad = (dm[[1, "Ga"]] - dm[[2, "Ga"]])/(2 h);
    <|"d_step" -> h, "Gd" -> gd, "Gad" -> gad,
      "Gd_second" -> (dm[[1, "G"]] + dm[[2, "G"]] - 2 z["G"])/h^2,
      "Ga_dplus" -> dm[[1, "Ga"]], "Ga_dminus" -> dm[[2, "Ga"]]|>, {h, hp}];
  krows = Table[
    <|"v" -> rows[[i, "v"]], "Gvv" -> rows[[i, "Gvv"]], "Gavv" -> rows[[i, "Gavv"]],
      "radial_d_step" -> gdrows[[-1, "d_step"]], "Gd" -> gdrows[[-1, "Gd"]],
      "Gad" -> gdrows[[-1, "Gad"]],
      "K_radially_corrected" -> rows[[i, "Gavv"]] - gdrows[[-1, "Gad"]] rows[[i, "Gvv"]]/gdrows[[-1, "Gd"]],
      "d2_from_Gvv" -> -rows[[i, "Gvv"]]/(2 gdrows[[-1, "Gd"]])|>, {i, Length[rows]}];
  krows = Map[Join[#, <|"a2_from_K" -> -#["K_radially_corrected"]/(2 z["Gaa"]),
      "radial_quartic_delta" -> #["K_radially_corrected"]^2/(8 z["Gaa"] #["Gd"] )|>] &, krows];
  krows = Map[Join[#, <|
      "K_with_response_GdGad" -> #["Gavv"] - direct["Gad_response_B_kplus_splus"] #["Gvv"]/direct["Gd_response_B_1_kplus"],
      "K_response_GadGvv" -> #["Gavv"] - direct["Gad_response_B_kplus_splus"] direct["Gvv_response_minus_B_kminus_kminus"]/direct["Gd_response_B_1_kplus"]|>] &, krows];
  primal = mo["baseline_solution"]["raw_integrals"];
  gdWard = 2./mo["baseline_solution"]["lambda_over_P"];
  tplusErr = Abs[z["Gaa"] - primal["Tplus"]]; tplusPass = TrueQ[tplusErr < 10.^-8];
  ga0Small = TrueQ[Abs[z["Ga"]] < 10.^-8];
  <|"n" -> nn, "c_stationary" -> mo["c"], "baseline" -> z,
    "source_log_scale" -> mo["log_scale"], "source_tau" -> mo["tau"], "finite_window" -> {0., mo["hi"]},
    "outer_nodes" -> Length[mo["data"]], "soft_order" -> ord, "root_bisections" -> 48,
    "angular_steps" -> rows, "radial_steps" -> gdrows, "radially_corrected" -> krows,
    "direct_complete_response" -> direct,
    "Faa0_vs_primal_Tplus" -> <|"frozen_Gaa0" -> z["Gaa"], "primal_Tplus" -> primal["Tplus"],
      "absolute_error" -> tplusErr, "tolerance" -> 10.^-8, "pass" -> tplusPass|>,
    "explicit_Gaa0_and_Ga0_assertions" -> <|"Gaa0_equals_Tplus_abs_lt_1e-8" -> tplusPass,
      "Ga0_abs_lt_1e-8" -> ga0Small, "Ga0" -> z["Ga"]|>,
    "Gd_vs_Ward_2_over_lambdaOverP" -> <|"Ward_value" -> gdWard,
      "fd_d001" -> gdrows[[1, "Gd"]], "fd_d0005" -> gdrows[[2, "Gd"]],
      "response_direct" -> direct["Gd_response_B_1_kplus"],
      "absolute_error_d0005" -> Abs[gdrows[[2, "Gd"]] - gdWard],
      "response_minus_Ward" -> direct["Gd_response_B_1_kplus"] - gdWard,
      "response_Gad_minus_fd" -> direct["Gad_response_B_kplus_splus"] - gdrows[[-1, "Gad"]]|>,
    "baseline_mean_residual" -> mo["baseline_solution"]["mean_residual"]|>
];

If[!TrueQ[$compute58LibraryMode],
start58 = AbsoluteTime[]; models58 = {};
base58[runStatus_] := <|
  "status" -> runStatus, "wolfram" -> $Version,
  "processor_count_reported" -> $ProcessorCount,
  "kernel_concurrency_requested" -> 1,
  "started" -> DateString[DateList[start58], {"ISODate", " ", "Time"}],
  "elapsed_seconds" -> N[AbsoluteTime[] - start58],
  "definition" -> "G(d,v,a)=E_{c+a,c+a}[H0(T1+d+v,T2+d-v)], with H0/lambda/roots/tau fixed from exactly one source solve52 baseline; receiver d,v,a evaluations never reoptimize the source",
  "correction" -> "K=Gavv-Gad*Gvv/Gd; a2=-K/(2Gaa); radial_quartic_delta=K^2/(8 Gaa Gd)",
  "models" -> models58,
  "machine_warnings" -> {"Finite square and machine quadrature only; no interval enclosure.", "Inherited expSafe52 drops individual log-weight contributions below -700.", "Positive Gamma upper-tail diagnostic uses GammaRegularized[...,x,Infinity], never 1-CDF.", "Angular/radial derivatives are centered finite differences; complete-response values are a separate quadrature cross-check.", "No rigorous sign claim for K and no C_N-gain claim."},
  "scope" -> "Frozen-source local micro-model diagnostic, not a proof about the general strong maximal constant."
|>;
save58[runStatus_] := Module[{x = base58[runStatus]}, If[Export[out58, x, "RawJSON"] === $Failed, Print["Export failed"]; Exit[3]]; If[!AssociationQ[Quiet@Check[Import[out58, "RawJSON"], $Failed]], Print["JSON validation failed"]; Exit[4]]];
run58[nn_, no_, ord_, tag_] := Module[{t0 = AbsoluteTime[], rr, pass58},
  Print[<|"starting" -> tag, "N" -> nn, "outer_nodes" -> no, "soft_order" -> ord|>];
  rr = response57[nn, no, ord];
  rr = Join[rr, <|"tag" -> tag, "model_seconds" -> N[AbsoluteTime[] - t0], "completed" -> DateString[DateList[], {"ISODate", " ", "Time"}]|>];
  AppendTo[models58, rr];
  pass58 = TrueQ[rr["explicit_Gaa0_and_Ga0_assertions", "Gaa0_equals_Tplus_abs_lt_1e-8"]] &&
    TrueQ[rr["explicit_Gaa0_and_Ga0_assertions", "Ga0_abs_lt_1e-8"]] &&
    TrueQ[Abs[rr["baseline_mean_residual"]] < 10.^-9];
  If[!pass58, save58["failed validation after " <> tag]; Print[<|"validation_failed" -> tag|>]; Exit[5]];
  save58["partial; completed " <> tag];
  Print[<|"completed" -> tag, "seconds" -> rr["model_seconds"], "K_v0125" -> rr["radially_corrected"][[-1, "K_radially_corrected"]]|>]
];

(* The fine replay targets N=1024, as the preselected highest-dimensional refinement. *)
save58["partial; initialized"];
run58[512, 256, 32, "N512_coarse_256x32"];
run58[1024, 256, 32, "N1024_coarse_256x32"];
run58[1024, 512, 64, "N1024_selected_512x64"];
save58["complete; all requested models exported"];
Print[<|"status" -> "complete", "seconds" -> N[AbsoluteTime[] - start58], "output" -> out58|>];
];
