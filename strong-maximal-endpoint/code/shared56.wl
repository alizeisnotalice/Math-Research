(* Exact finite layer-cake certificate for the two-slab common-action example. *)
Clear[t];
payoff = {{1/2, 0}, {1/4, 1/2}}; (* endpoint matrix for the common-average objective *)
adaptive = Mean[Max /@ payoff];
shared = Max[Mean /@ Transpose[payoff]];
(* For the actual continuum of intervals, t=3 covers both fibers up to level 1/3. *)
layerGap[t_] := Piecewise[{{0, 0 < t <= 1/3}, {1/2, 1/3 < t < 1/2}, {0, t >= 1/2}}];
Print["Wolfram version: ", $Version];
Print["adaptive, shared, regret = ", {adaptive, shared, adaptive - shared}];
Print["integrated layer gap = ", Integrate[layerGap[t], {t, 0, 1/2}]];
Print["regret minus layer certificate = ", FullSimplify[adaptive - shared - Integrate[layerGap[t], {t, 0, 1/2}]]];
g1[r_] := Piecewise[{{1, 0 < r <= 1}, {1/r, r > 1}}];
g2[r_] := Piecewise[{{1, 0 < r <= 1}, {3 - 2/r, 1 < r <= 2}, {4/r, r > 2}}];
Print["height-cover counterexample m1,m2,S,A = ",
 {MaxValue[{g1[r], r > 0}, r], MaxValue[{g2[r], r > 0}, r],
  MaxValue[{(g1[r] + g2[r])/2, r > 0}, r], (1 + 2)/2}];

actualVals={MaxValue[{g1[r],r>0},r],MaxValue[{g2[r],r>0},r],MaxValue[{(g1[r]+g2[r])/2,r>0},r]};
mask[z_]:=Piecewise[{{1,t<z}},0];
coverCheck=TrueQ[Resolve[ForAll[t,t>0,Max[2mask[1],mask[1/2]+mask[2]]==mask[1]+mask[2]],Reals]];
checks=<|"old_slab_regret"->TrueQ[adaptive-shared==1/8],"old_slab_continuous_layer_integral"->TrueQ[Integrate[layerGap[t],{t,0,1/2}]==1/12],"new_actual_fiber_and_shared_max"->TrueQ[actualVals=={1,2,5/4}],"new_two_actions_cover_every_height"->coverCheck,"new_positive_regret"->TrueQ[3/2-actualVals[[3]]==1/4]|>;
If[!And@@Values[checks],Print[checks];Exit[2]];
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/shared56.json"];
Export[out,<|"wolfram"->$Version,"checks"->checks,"all_pass"->True,"scope"->"exact fixed-outer-interval counterexample; not a halo statement","actual_Phi2_budget"->ToString[2+3(1+Log[3]),InputForm]|>,"RawJSON"];
