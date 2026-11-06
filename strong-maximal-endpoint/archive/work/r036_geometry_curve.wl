(* r036: two genuinely deformed two-block receiver curves with joint source
   reweighting.  Direct split quadrature only; floating evidence, no FOSD. *)
root=DirectoryName[$InputFileName];args=$ScriptCommandLine;repo=FileNameJoin[{DirectoryName[root],"Math-Research-strong-maximal","strong-maximal-endpoint"}];
seedFile=If[Length[args]>=2,args[[2]],FileNameJoin[{repo,"rounds","r030","receiver_seed.json"}]];out=If[Length[args]>=3,args[[3]],FileNameJoin[{root,"r036_geometry_curve.json"}]];
$r032LibraryMode=True;Get[FileNameJoin[{repo,"code","receiver_continuous32.wl"}]];$r032LibraryMode=False;
dBase=dd;wBase=w0;srcScale=1.0276924765700282;recC=c+.1075;
makeU[dnew_]:=Join[{{0.,0.}},Transpose[{Rest[dnew]+Rest[vv],Rest[dnew]-Rest[vv]}],Transpose[{Rest[dnew]-Rest[vv],Rest[dnew]+Rest[vv]}]];
coneVolume[uu_]:=Module[{s=SortBy[uu,First],xx,best,tot=0.,i},xx=DeleteDuplicates[s[[All,1]]];Do[best=Min[Select[s,#[[1]]<=xx[[i]]&][[All,2]]];tot+=(Exp[-xx[[i]]]-If[i<Length[xx],Exp[-xx[[i+1]]],0.])Exp[-best],{i,Length[xx]}];tot];
(* First-order mass-balanced transfer between inner (v=2..5) and outer
   (v=9.5,10,10.5,16) source layers; it is not exact mass preservation. *)
inner={5,6,7,8,9,10,11};outer={19,20,21,33};si=Total[wBase[[inner]]];so=Total[wBase[[outer]]];
reweight[a_]:=Module[{z=wBase},z[[inner]]*=Exp[a];z[[outer]]*=Exp[-a si/so];z];
setGeometry[dnew_,ww_]:=(u=makeU[dnew];lo=Min[u];hi=Max[u]+b mz+16. Sqrt[b vz];w0=ww;);
acceptedAnchors[ids_]:=Flatten[Table[If[j==1,{u[[1]]},{u[[j]],u[[j+ng-1]]}],{j,ids}],1];
compact[r_]:=Module[{av=r["averages_by_group"],ids,an,vv0},ids=Flatten@Position[av,_?(#>1.00001&)];an=acceptedAnchors[ids];vv0=If[Length[an]==0,0.,coneVolume[an]];<|"scaled_cost"->r["scaled_cost"],"actual_logJ"->r["actual_logJ"],"min_average"->r["min_average"],"min_group"->First@Ordering[av,1],"accepted_groups_margin_1e5"->ids,"accepted_anchor_count"->Length[an],"accepted_cone_volume"->vv0,"accepted_cone_log_score"->If[vv0>0,Log[vv0]-r["actual_logJ"],-Infinity],"averages_groups_22_to_26"->av[[22;;26]],"seconds"->r["seconds"]|>];
run[name_,dnew_,a_,no_:128,ord_:32]:=Module[{r},setGeometry[dnew,reweight[a]];r=evaluateReceiverC[ConstantArray[recC,ng],srcScale,no,ord];<|"geometry"->name,"alpha"->a,"d_min"->Min[dnew],"d_max"->Max[dnew],"result"->compact[r]|>];
v=vv;dOuter=dBase+.15 (v/16.)^2;dMiddle=dBase-.10 v (16.-v)/64.;
base=run["undeformed",dBase,0.];trials=Join[Table[run["outer_up",dOuter,a],{a,{-.015,0.,.015}}],Table[run["middle_in",dMiddle,a],{a,{-.015,0.,.015}}]];
best=First@Ordering[trials[[All,"result","min_average"]],-1];winner=trials[[best]];
winner32=run[winner["geometry"],If[winner["geometry"]==="outer_up",dOuter,dMiddle],winner["alpha"],256,32];
winner64=run[winner["geometry"],If[winner["geometry"]==="outer_up",dOuter,dMiddle],winner["alpha"],256,64];
vBase=base["result","accepted_cone_volume"];vWinner=winner64["result","accepted_cone_volume"];
scoreBase=base["result","accepted_cone_log_score"];scoreWinner=winner64["result","accepted_cone_log_score"];scalarScore=71.71499208794934;
(* restore no shared-state surprise for callers *)
setGeometry[dBase,wBase];
res=<|"status"->"r036 deformed receiver-curve plus joint source reweighting; floating split-quadrature evidence only","n"->n,"source_scale"->srcScale,"global_receiver_c"->recC,
 "curve_formulas"-><|"outer_up"->"d(v)+0.15(v/16)^2","middle_in"->"d(v)-0.10 v(16-v)/64"|>,"reweight"->"inner active groups 5..11 multiply exp(alpha); outer active groups 19,20,21,33 multiply exp(-alpha sum_inner/sum_outer), a first-order total-mass balance only",
 "baseline_N128_soft32"->base,"trials_N128_soft32"->trials,"best_trial"->winner,"best_trial_N256_soft32"->winner32,"best_trial_N256_soft64"->winner64,
 "finite_pareto_cone_volume"-><|"baseline"->vBase,"best_trial"->vWinner,"baseline_log_score"->scoreBase,"best_trial_log_score"->scoreWinner,"best_minus_baseline"->scoreWinner-scoreBase,"scalar_log_score"->scalarScore,"best_minus_scalar"->scoreWinner-scalarScore|>,
 "same_n_scalar_reference_log_score"->scalarScore,"same_n_scalar_reference_scope"->"r026 n=64 same-c, height-2 numerical scalar benchmark",
 "interpretation"->"Each trial rebuilds source H from changed curve u and changed source weights before integrating true Phi_n cost. The finite Pareto cone includes only groups with numerical average above 1+1e-5 (and both mirror anchors); no high-dimensional lift is inferred."|>;
Export[out,res,"RawJSON"];Print[<|"base_min"->base["result","min_average"],"best"->winner|>];
