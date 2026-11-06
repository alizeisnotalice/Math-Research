(* r057: frozen-source, radial-compensated receiver response diagnostic.
   Floating finite-window quadrature only; source profile, roots, tau and
   calibrated lambda remain fixed while (d,v,a) vary. *)
root57 = DirectoryName[$InputFileName]; args57 = $ScriptCommandLine;
out57 = If[Length[args57] >= 2, args57[[2]], FileNameJoin[{$TemporaryDirectory, "frozen57.json"}]];
env57=Quiet@Check[Environment["STRONG_MAXIMAL_REPO"],$Failed];
repo57=If[StringQ[env57]&&env57!="",env57,DirectoryName[root57]];
If[!FileExistsQ[FileNameJoin[{repo57,"code","primal52.wl"}]],Print["Set STRONG_MAXIMAL_REPO to project root"];Exit[2]];
Block[{$ScriptCommandLine = {}, $primal52LibraryMode = True},
  Get[FileNameJoin[{repo57, "code", "primal52.wl"}]]
];
seed57 = Import[FileNameJoin[{repo57, "rounds", "r052", "primal52.json"}], "RawJSON"];
seedRow57[nn_] := SelectFirst[seed57["rows"], #["n"] == nn &];

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

response57[nn_, no_Integer, ord_Integer] := Module[
  {mo, z, hp = {.01, .005}, vs = {.5, .25, .125}, vm, dm, gdrows, gd, gad,
   rows, krows, primal, h, gdWard, tplusErr, tplusPass},
  mo = makeFrozen57[nn, no, ord]; z = evalFrozen57[mo, 0., 0., 0.];
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
  primal = mo["baseline_solution"]["raw_integrals"];
  gdWard = 2./mo["baseline_solution"]["lambda_over_P"];
  tplusErr = Abs[z["Gaa"] - primal["Tplus"]]; tplusPass = TrueQ[tplusErr < 10.^-8];
  <|"n" -> nn, "c_stationary" -> mo["c"], "baseline" -> z,
    "source_log_scale" -> mo["log_scale"], "source_tau" -> mo["tau"], "finite_window" -> {0., mo["hi"]},
    "outer_nodes" -> Length[mo["data"]], "soft_order" -> ord, "root_bisections" -> 48,
    "angular_steps" -> rows, "radial_steps" -> gdrows, "radially_corrected" -> krows,
    "Faa0_vs_primal_Tplus" -> <|"frozen_Gaa0" -> z["Gaa"], "primal_Tplus" -> primal["Tplus"],
      "absolute_error" -> tplusErr, "tolerance" -> 10.^-8, "pass" -> tplusPass|>,
    "Gd_vs_Ward_2_over_lambdaOverP" -> <|"Ward_value" -> gdWard,
      "fd_d001" -> gdrows[[1, "Gd"]], "fd_d0005" -> gdrows[[2, "Gd"]],
      "absolute_error_d0005" -> Abs[gdrows[[2, "Gd"]] - gdWard]|>,
    "baseline_mean_residual" -> mo["baseline_solution"]["mean_residual"],
    "omitted_receiver_probability_union_diagnostic" -> 2 Dot[weights52[mo["b"],mo["c"]],GammaRegularized[mo["b"],Clip[mo["hi"]-.51-mo["c"] Range[0,mo["b"]],{0.,Infinity}],Infinity]],
    "tail_scope" -> "Positive Gamma upper-tail diagnostic covers |d|<=.01, |v|<=.5, receiver a offset zero; not a certified derivative-tail enclosure"|>
];

start57 = AbsoluteTime[];
coarse57 = response57[64, 256, 32]; coarse12857 = response57[128, 256, 32]; coarse25657 = response57[256, 256, 32];
fine6457 = response57[64, 512, 64]; fine25657 = response57[256, 512, 64];
res57 = <|
  "status" -> "floating frozen-source radial-compensation diagnostic; not an interval certificate",
  "wolfram" -> $Version, "seconds" -> N[AbsoluteTime[] - start57],
  "definition" -> "G(d,v,a)=E_{c+a,c+a}[H0(T1+d+v,T2+d-v)], with H0/lambda/roots fixed from solve52 baseline",
  "correction" -> "K=Gavv-Gad*Gvv/Gd, not raw Gavv; a2=-K/(2Gaa), radial_quartic_delta=K^2/(8 Gaa Gd)",
  "coarse_256x32" -> {coarse57, coarse12857, coarse25657},
  "validation_512x64" -> {fine6457, fine25657},
  "scope" -> "Finite square and inherited expSafe52 truncation. Analytic receiver c derivatives use cdfJets52 and jetBuild52; radial/v derivatives are centered finite differences. No optimizer or certificate claim."
|>;
If[Export[out57,res57,"RawJSON"]===$Failed,Print["Export failed"];Exit[3]];If[!AssociationQ[Import[out57,"RawJSON"]],Exit[4]];Print[<|"status"->res57["status"],"seconds"->res57["seconds"],"output"->out57|>];
