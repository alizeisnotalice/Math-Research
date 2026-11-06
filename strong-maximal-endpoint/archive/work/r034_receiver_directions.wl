(* r034: direct mirrored two-coordinate receiver probes about r033's global c.
   Uses the repository r032 split integrator as a library; no FOSD assumption. *)
root34 = DirectoryName[$InputFileName]; args34 = $ScriptCommandLine;
repo34 = FileNameJoin[{DirectoryName[root34], "Math-Research-strong-maximal", "strong-maximal-endpoint"}];
seed34 = If[Length[args34] >= 2, args34[[2]], FileNameJoin[{repo34, "rounds", "r030", "receiver_seed.json"}]];
out34 = If[Length[args34] >= 3, args34[[3]], FileNameJoin[{root34, "r034_receiver_directions.json"}]];
$r032LibraryMode = True; Get[FileNameJoin[{repo34, "code", "receiver_continuous32.wl"}]]; $r032LibraryMode = False;

(* For j>1: u[j] has coordinate parameters (cA,cB), and its reflected point
   u[j+ng-1] has (cB,cA).  This is the correct mirror-swapped convention. *)
prob34[t1_?NumericQ,t2_?NumericQ,ca_List,cb_List] := Module[{ans=ConstantArray[0.,ng],z,q1,q2},
 ans[[1]]=polyC[t1,ca[[1]]] polyC[t2,cb[[1]]] Exp[-t1-t2];
 Do[z=u[[j]];q1=polyC[t1-z[[1]],ca[[j]]];q2=polyC[t2-z[[2]],cb[[j]]];ans[[j]]+=q1 q2 Exp[z[[1]]+z[[2]]-t1-t2];
   z=u[[j+ng-1]];q1=polyC[t1-z[[1]],cb[[j]]];q2=polyC[t2-z[[2]],ca[[j]]];ans[[j]]+=q1 q2 Exp[z[[1]]+z[[2]]-t1-t2],{j,2,ng}];ans];
constant34[t1_?NumericQ,l_?NumericQ,r_?NumericQ,h_?NumericQ,ca_List,cb_List] := Module[{cost=0.,av=ConstantArray[0.,ng],z,fac},
 If[r<=l,Return[{cost,av}]];cost=Exp[tau-t1](Exp[-l]-Exp[-r])phi[h];
 fac=polyC[t1,ca[[1]]] Exp[-t1];av[[1]]=fac(gbcdfC[r,cb[[1]]]-gbcdfC[l,cb[[1]]]);
 Do[z=u[[j]];fac=polyC[t1-z[[1]],ca[[j]]]Exp[z[[1]]-t1];av[[j]]+=fac(gbcdfC[r-z[[2]],cb[[j]]]-gbcdfC[l-z[[2]],cb[[j]]]);
   z=u[[j+ng-1]];fac=polyC[t1-z[[1]],cb[[j]]]Exp[z[[1]]-t1];av[[j]]+=fac(gbcdfC[r-z[[2]],ca[[j]]]-gbcdfC[l-z[[2]],ca[[j]]]),{j,2,ng}];{cost,av}];
soft34[t1_?NumericQ,l_?NumericQ,r_?NumericQ,amp_,ca_List,cb_List,ord_Integer] := Module[{x,wt,cost=0.,av=ConstantArray[0.,ng],hh},
 If[r<=l,Return[{cost,av}]];{x,wt}=quad[ord,l,r];
 Do[hh=profile[likeAmp[x[[i]],amp]];cost+=wt[[i]]Exp[tau-t1-x[[i]]]phi[hh];av+=wt[[i]]hh prob34[t1,x[[i]],ca,cb],{i,Length[x]}];{cost,av}];
run34[ca_List,cb_List,no_Integer,ord_Integer] := Module[{ww=srcScale w0,tx,tw,cost=0.,av=ConstantArray[0.,ng],amp,roots,z,a,start=AbsoluteTime[]},
 {tx,tw}=outerQuad[no];Do[amp=sourceAmp[tx[[i]],ww];roots=boundary[#,amp]&/@{1.,N[n],phiPrime};
  {z,a}=constant34[tx[[i]],roots[[1]],roots[[2]],1.,ca,cb];cost+=tw[[i]]z;av+=tw[[i]]a;
  {z,a}=constant34[tx[[i]],roots[[3]],hi,cap,ca,cb];cost+=tw[[i]]z;av+=tw[[i]]cap a;
  {z,a}=soft34[tx[[i]],roots[[2]],roots[[3]],amp,ca,cb,ord];cost+=tw[[i]]z;av+=tw[[i]]a,{i,Length[tx]}];
 <|"outer_nodes"->Length[tx],"soft_inner_order"->ord,"scaled_cost"->cost,"actual_logJ"->Log[cost]-tau,
  "averages_by_group"->av/coeff,"min_average"->Min[av/coeff],"min_group"->First@Ordering[av/coeff,1],"seconds"->AbsoluteTime[]-start|>];

srcScale=1.0276924765700282; globalC=c+.10; jstar=24; baseA=ConstantArray[globalC,ng];baseB=ConstantArray[globalC,ng];
mk[da_,db_,no_:128,ord_:32] := Module[{aa=baseA,bb=baseB,rr},aa[[jstar]]+=da;bb[[jstar]]+=db;rr=run34[aa,bb,no,ord];<|"delta_c1"->da,"delta_c2"->db,"c1"->aa[[jstar]],"c2"->bb[[jstar]],"result"->rr|>];
baseRun=mk[0.,0.]; steps={.005,.01,.0125,.025,.05}; probes=Join[Flatten[Table[{mk[h,0.],mk[-h,0.],mk[0.,h],mk[0.,-h]},{h,steps}],1]];
summary[x_] := <|"d1"->x["delta_c1"],"d2"->x["delta_c2"],"min"->x["result","min_average"],"min_group"->x["result","min_group"],"group24_average"->x["result","averages_by_group"][[jstar]]|>;
bestIndex=First@Ordering[(#["result","min_average"]& /@ probes),-1];best=probes[[bestIndex]];
(* Only a genuinely better same-grid minimum earns a 256 validation. *)
validation=If[best["result","min_average"]>baseRun["result","min_average"]+10^-10,
 <|"N256_soft32"->mk[best["delta_c1"],best["delta_c2"],256,32],"N256_soft64"->mk[best["delta_c1"],best["delta_c2"],256,64]|>,
 <|"not_run_reason"->"no N128 perturbation exceeds same-grid global baseline"|>];
grad25=<|"d_c1_central_h005"->(probes[[1,"result","averages_by_group",jstar]]-probes[[2,"result","averages_by_group",jstar]])/.01,
 "d_c2_central_h005"->(probes[[3,"result","averages_by_group",jstar]]-probes[[4,"result","averages_by_group",jstar]])/.01|>;
out=<|"status"->"r034 direct mirrored-coordinate perturbations; floating numerical evidence only", "source_scale"->srcScale,"global_receiver_c"->globalC,
 "group"->jstar,"mapping"->"u[j]:(c1,c2); u[j+32]:(c2,c1)","baseline_N128_soft32"->baseRun,
 "probes_N128_soft32"->probes,"compact_probe_summary"->(summary /@ probes),"central_gradient_h005_group24"->grad25,
 "best_N128"->best,"validation"->validation,"interpretation"->"Direct quadrature only. The minimum across groups can be nonsmooth; no stochastic-order claim is made."|>;
Export[out34,out,"RawJSON"];Print[<|"baseline"->baseRun["result","min_average"],"best"->best["result","min_average"],"best_d"->{best["delta_c1"],best["delta_c2"]}|>];
