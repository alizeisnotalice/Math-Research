(* r036: two genuinely deformed two-block receiver curves with joint source
   reweighting.  Direct split quadrature only; floating evidence, no FOSD. *)
root=DirectoryName[$InputFileName];args=$ScriptCommandLine;repo=DirectoryName[root];
seedFile=If[Length[args]>=2,args[[2]],FileNameJoin[{repo,"rounds","r030","receiver_seed.json"}]];outFile36=If[Length[args]>=3,args[[3]],FileNameJoin[{$TemporaryDirectory,"geometry_probe36.json"}]];
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
bestScoreTrial=trials[[First@Ordering[trials[[All,"result","accepted_cone_log_score"]],-1]]];
winner32=run[winner["geometry"],If[winner["geometry"]==="outer_up",dOuter,dMiddle],winner["alpha"],256,32];
winner64=run[winner["geometry"],If[winner["geometry"]==="outer_up",dOuter,dMiddle],winner["alpha"],256,64];
vBase=base["result","accepted_cone_volume"];vWinner=winner64["result","accepted_cone_volume"];
scoreBase=base["result","accepted_cone_log_score"];scoreWinner=winner64["result","accepted_cone_log_score"];scalarScore=71.71499208794934;
(* restore no shared-state surprise for callers *)
setGeometry[dBase,wBase];
res=<|"status"->"r036 deformed receiver-curve plus joint source reweighting; floating split-quadrature evidence only","n"->n,"source_scale"->srcScale,"global_receiver_c"->recC,
 "curve_formulas"-><|"outer_up"->"d(v)+0.15(v/16)^2","middle_in"->"d(v)-0.10 v(16-v)/64"|>,"reweight"->"inner active groups 5..11 multiply exp(alpha); outer active groups 19,20,21,33 multiply exp(-alpha sum_inner/sum_outer), a first-order total-mass balance only",
 "baseline_N128_soft32"->base,"trials_N128_soft32"->trials,"highest_min_average_trial"->winner,"best_score_trial_N128_soft32"->bestScoreTrial,"best_score_minus_baseline"->bestScoreTrial["result","accepted_cone_log_score"]-scoreBase,"highest_min_trial_N256_soft32"->winner32,"highest_min_trial_N256_soft64"->winner64,
 "finite_pareto_cone_volume"-><|"baseline"->vBase,"highest_min_trial"->vWinner,"baseline_log_score"->scoreBase,"highest_min_trial_log_score"->scoreWinner,"highest_min_trial_minus_baseline"->scoreWinner-scoreBase,"scalar_log_score"->scalarScore,"highest_min_trial_minus_scalar"->scoreWinner-scalarScore|>,
 "same_n_scalar_reference_log_score"->scalarScore,"same_n_scalar_reference_scope"->"r026 n=64 height-2 scalar numerical reference at seed c0; not the full Q64 optimum",
 "interpretation"->"Each trial rebuilds source H from changed curve u and changed source weights before integrating true Phi_n cost. The finite Pareto cone includes only groups with numerical average above 1+1e-5 (and both mirror anchors); no high-dimensional lift is inferred."|>;
geometryResult36=res;

srcScale=1.0276924765700282;recC=c+.1075;dnew=dd-.10 vv(16.-vv)/64.;wBase=w0;inner={5,6,7,8,9,10,11};outer={19,20,21,33};si=Total[wBase[[inner]]];so=Total[wBase[[outer]]];
wSrc=Module[{z=wBase,a=.015},z[[inner]]*=Exp[a];z[[outer]]*=Exp[-a si/so];z];u=Join[{{0.,0.}},Transpose[{Rest[dnew]+Rest[vv],Rest[dnew]-Rest[vv]}],Transpose[{Rest[dnew]-Rest[vv],Rest[dnew]+Rest[vv]}]];lo=Min[u];hi=Max[u]+b mz+16.Sqrt[b vz];w0=wSrc;
cv[uu_]:=Module[{s=SortBy[uu,First],xx,best,tot=0.,i},xx=DeleteDuplicates[s[[All,1]]];Do[best=Min[Select[s,#[[1]]<=xx[[i]]&][[All,2]]];tot+=(Exp[-xx[[i]]]-If[i<Length[xx],Exp[-xx[[i+1]]],0.])Exp[-best],{i,Length[xx]}];tot];
pv[t1_,t2_,aa_,bb_]:=Module[{ans=ConstantArray[0.,ng],z,q1,q2},ans[[1]]=polyC[t1,aa[[1]]]polyC[t2,bb[[1]]]Exp[-t1-t2];Do[z=u[[j]];q1=polyC[t1-z[[1]],aa[[j]]];q2=polyC[t2-z[[2]],bb[[j]]];ans[[j]]+=q1 q2 Exp[z[[1]]+z[[2]]-t1-t2];z=u[[j+ng-1]];q1=polyC[t1-z[[1]],bb[[j]]];q2=polyC[t2-z[[2]],aa[[j]]];ans[[j]]+=q1 q2 Exp[z[[1]]+z[[2]]-t1-t2],{j,2,ng}];ans];
cp[t1_,l_,r_,h_,aa_,bb_]:=Module[{z,fac,cost=0.,av=ConstantArray[0.,ng]},If[r<=l,Return[{cost,av}]];cost=Exp[tau-t1](Exp[-l]-Exp[-r])phi[h];fac=polyC[t1,aa[[1]]]Exp[-t1];av[[1]]=fac(gbcdfC[r,bb[[1]]]-gbcdfC[l,bb[[1]]]);Do[z=u[[j]];fac=polyC[t1-z[[1]],aa[[j]]]Exp[z[[1]]-t1];av[[j]]+=fac(gbcdfC[r-z[[2]],bb[[j]]]-gbcdfC[l-z[[2]],bb[[j]]]);z=u[[j+ng-1]];fac=polyC[t1-z[[1]],bb[[j]]]Exp[z[[1]]-t1];av[[j]]+=fac(gbcdfC[r-z[[2]],aa[[j]]]-gbcdfC[l-z[[2]],aa[[j]]]),{j,2,ng}];{cost,av}];
sp[t1_,l_,r_,amp_,aa_,bb_]:=Module[{x,wt,h,cost=0.,av=ConstantArray[0.,ng]},If[r<=l,Return[{cost,av}]];{x,wt}=quad[32,l,r];Do[h=profile[likeAmp[x[[i]],amp]];cost+=wt[[i]]Exp[tau-t1-x[[i]]]phi[h];av+=wt[[i]]h pv[t1,x[[i]],aa,bb],{i,Length[x]}];{cost,av}];
run[aa_,bb_]:=Module[{ww=srcScale w0,tx,tw,amp,roots,z,a,cost=0.,av=ConstantArray[0.,ng],ids,an,vol},{tx,tw}=outerQuad[128];Do[amp=sourceAmp[tx[[i]],ww];roots=boundary[#,amp]&/@{1.,N[n],phiPrime};{z,a}=cp[tx[[i]],roots[[1]],roots[[2]],1.,aa,bb];cost+=tw[[i]]z;av+=tw[[i]]a;{z,a}=cp[tx[[i]],roots[[3]],hi,cap,aa,bb];cost+=tw[[i]]z;av+=tw[[i]]cap a;{z,a}=sp[tx[[i]],roots[[2]],roots[[3]],amp,aa,bb];cost+=tw[[i]]z;av+=tw[[i]]a,{i,Length[tx]}];ids=Flatten@Position[av/coeff,_?(#>1.00001&)];an=Flatten[Table[If[j==1,{u[[1]]},{u[[j]],u[[j+ng-1]]}],{j,ids}],1];vol=If[Length[an]>0,cv[an],0.];<|"scaled_cost"->cost,"actual_logJ"->Log[cost]-tau,"min_average"->Min[av/coeff],"accepted_groups_margin_1e5"->ids,"accepted_cone_volume"->vol,"accepted_cone_log_score"->If[vol>0,Log[vol]-(Log[cost]-tau),-Infinity]|>];
q=vv/16.;baseA=ConstantArray[recC,ng];baseB=ConstantArray[recC,ng];cases={<|"name"->"uniform_recC","c1"->baseA,"c2"->baseB|>,<|"name"->"mirror_kappa_plus","c1"->recC+.1q,"c2"->recC-.1q|>,<|"name"->"mirror_kappa_minus","c1"->recC-.1q,"c2"->recC+.1q|>,<|"name"->"common_quadratic","c1"->recC+.1q^2,"c2"->recC+.1q^2|>};
rows=(<|"name"->#["name"],"result"->run[#["c1"],#["c2"]]|>& /@ cases);res=<|"status"->"r036 receiver-dependent directions on fixed middle-in alpha=.015 source; floating direct split quadrature only","source_geometry"->"d(v)-.10v(16-v)/64","source_reweight_alpha"->.015,"mapping"->"u[j]:(c1,c2), mirror u[j+32]:(c2,c1)","threshold"->"accepted avg>1+1e-5","rows_N128_soft32"->rows,"interpretation"->"All rows share the same source profile; no FOSD or scalar-Q conclusion is claimed."|>;
receiverResult36=res;
Export[outFile36,<|"wolfram"->$Version,"input_dimension"->n,"block_size"->b,"status"->"continuous split quadrature, finite accepted cone volume; not interval certified","geometry_joint"->geometryResult36,"receiver_dependent"->receiverResult36|>,"RawJSON"];
Print[<|"output"->outFile36,"joint_best_minus_baseline"->geometryResult36["finite_pareto_cone_volume","highest_min_trial_minus_baseline"],"dependent_scores"->receiverResult36["rows_N128_soft32"][[All,"result","accepted_cone_log_score"]]|>];
