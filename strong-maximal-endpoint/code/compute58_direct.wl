(* Direct K replay from the immutable, completed N1024 512x64 frozen-source
   record.  It rebuilds only quadrature nodes and the already-recorded source
   level sets; it never calls solve52 or recalibrates lambda. *)
root58f = DirectoryName[$InputFileName]; out58f = If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{$TemporaryDirectory,"compute58.json"}]];
$compute58LibraryMode = True; Get[FileNameJoin[{root58f, "compute58.wl"}]]; $compute58LibraryMode = False;
start58f = AbsoluteTime[]; If[!FileExistsQ[out58f],Print["Run compute58.wl first; missing input JSON: "<>out58f];Exit[2]]; old58f = Import[out58f,"RawJSON"];If[!AssociationQ[old58f],Exit[2]];
rec58f = SelectFirst[old58f["models"], #["tag"] == "N1024_selected_512x64" &,Missing["record"]];If[MissingQ[rec58f],Print["Selected source record missing"];Exit[2]];
n58f = rec58f["n"]; c58f = rec58f["c_stationary"]; hi58f = rec58f["finite_window"][[2]]; no58f = rec58f["outer_nodes"]; ord58f = rec58f["soft_order"];
ctx58f = context52[n58f, c58f, c58f, no58f, ord58f, hi58f];
tauGap58f = Abs[ctx58f["tau"] - rec58f["source_tau"]];
If[!TrueQ[tauGap58f < 10.^-8], Print[<|"tau_record_mismatch" -> tauGap58f|>]; Exit[5]];
ls58f = rec58f["source_log_scale"]; js58f = ctx58f["j1"]; dat58f = Table[
  p58f = js58f[ctx58f["x"][[i]]]; la58f = ls58f - rec58f["source_tau"] + p58f[[5]];
  bd58f[lev_] := Module[{lo58f = 0., up58f = hi58f, md58f}, If[la58f + js58f[hi58f][[5]] < Log[lev], Return[hi58f]]; Do[md58f = (lo58f + up58f)/2.; If[la58f + js58f[md58f][[5]] < Log[lev], lo58f = md58f, up58f = md58f], {48}]; (lo58f + up58f)/2.];
  <|"x" -> ctx58f["x"][[i]], "w" -> ctx58f["w"][[i]], "logAmp" -> la58f, "roots" -> (bd58f /@ {1., N[n58f], ctx58f["aCap"]})|>,
  {i, Length[ctx58f["x"]]}];
mo58f = <|"n" -> n58f, "b" -> n58f/2, "c" -> c58f, "hi" -> hi58f, "source_j" -> js58f, "data" -> dat58f, "inner_rule" -> ctx58f["inner_rule"], "cap" -> ctx58f["cap"]|>;
st58f = sTCheck58[mo58f];
If[!TrueQ[st58f["absolute_difference"] < 10.^-8], Print[<|"s_t_validation_failed" -> st58f|>]; Exit[6]];
dk58f = directK58[mo58f]; fd58f = rec58f["radially_corrected"][[ -1 ]];
check58f = <|"record_tau_minus_rebuilt_tau" -> rec58f["source_tau"] - ctx58f["tau"], "FD_Gavv_v0125" -> fd58f["Gavv"], "direct_Gavv" -> dk58f["Gavv_direct_minus_B_kminus_stminus_plus_splus_kminus"], "direct_minus_FD_Gavv" -> dk58f["Gavv_direct_minus_B_kminus_stminus_plus_splus_kminus"] - fd58f["Gavv"], "FD_K_v0125" -> fd58f["K_radially_corrected"], "direct_K" -> dk58f["K_direct_response"], "direct_single_form_minus_direct_K" -> dk58f["K_direct_single_form"] - dk58f["K_direct_response"]|>;
tail58f[rr58f_] := Module[{bb58f = rr58f["n"]/2, cc58f = rr58f["c_stationary"], hh58f = rr58f["finite_window"][[2]], val58f},
  val58f = N[2. Dot[weights52[bb58f, cc58f], GammaRegularized[bb58f, Clip[hh58f - .51 - cc58f Range[0, bb58f], {0., Infinity}], Infinity]]];
  <|"tag" -> rr58f["tag"], "positive_Gamma_upper_tail" -> val58f, "scope" -> "covers |d|<=.01, |v|<=.5 and receiver a offset zero; it is not a certified derivative-tail enclosure"|>
];
tails58f = tail58f /@ old58f["models"];
new58f = Join[old58f, <|"positive_Gamma_upper_tail_diagnostics" -> tails58f, "direct_response_K_selected_N1024_512x64" -> <|"status" -> "floating direct complete-response cross-check; no interval sign conclusion", "source" -> "reconstructed from the completed N1024 512x64 record's c, log_scale, tau, fixed window and quadrature; no solve52/lambda recalibration", "s_t_analytic_vs_centered_fd" -> st58f, "response_values" -> dk58f, "comparison_to_completed_receiver_FD" -> check58f, "seconds" -> N[AbsoluteTime[] - start58f]|>|>];
If[Export[out58f, new58f, "RawJSON"] === $Failed, Exit[3]]; If[!AssociationQ[Import[out58f, "RawJSON"]], Exit[4]];
Print[KeyTake[check58f, {"direct_Gavv", "direct_minus_FD_Gavv", "direct_K", "FD_K_v0125"}]];
