(* Portable r034 fixed-source receiver-c audit.  Floating quadrature only. *)
here=DirectoryName[$InputFileName];args=$ScriptCommandLine;
roots=DeleteDuplicates[{DirectoryName[here],FileNameJoin[{DirectoryName[here],"Math-Research-strong-maximal","strong-maximal-endpoint"}]}];
repo=FirstCase[roots,x_/;FileExistsQ[FileNameJoin[{x,"code","receiver_continuous32.wl"}]]:>x,Missing["repo"]];
If[MissingQ[repo],Print["Cannot locate strong-maximal-endpoint repository relative to ",here];Exit[1]];
seedFile34=If[Length[args]>=2,args[[2]],FileNameJoin[{repo,"rounds","r030","receiver_seed.json"}]];
outFile34=If[Length[args]>=3,args[[3]],FileNameJoin[{$TemporaryDirectory,"r034_receiver_final.json"}]];
$r032LibraryMode=True;Get[FileNameJoin[{repo,"code","receiver_continuous32.wl"}]];$r032LibraryMode=False;

(* Correct two-coordinate mirror convention: u[j] uses (c1,c2), while
   u[j+ng-1] uses (c2,c1).  No stochastic-order claim is used. *)
pv[t1_?NumericQ,t2_?NumericQ,aa_List,bb_List] := Module[{z,q1,q2,outv=ConstantArray[0.,ng]},
 outv[[1]]=polyC[t1,aa[[1]]]polyC[t2,bb[[1]]]Exp[-t1-t2];
 Do[z=u[[j]];q1=polyC[t1-z[[1]],aa[[j]]];q2=polyC[t2-z[[2]],bb[[j]]];outv[[j]]+=q1 q2 Exp[z[[1]]+z[[2]]-t1-t2];
  z=u[[j+ng-1]];q1=polyC[t1-z[[1]],bb[[j]]];q2=polyC[t2-z[[2]],aa[[j]]];outv[[j]]+=q1 q2 Exp[z[[1]]+z[[2]]-t1-t2],{j,2,ng}];outv];
cp[t1_?NumericQ,l_?NumericQ,r_?NumericQ,h_?NumericQ,aa_List,bb_List] := Module[{z,fac,cost=0.,av=ConstantArray[0.,ng]},
 If[r<=l,Return[{cost,av}]];cost=Exp[tau-t1](Exp[-l]-Exp[-r])phi[h];fac=polyC[t1,aa[[1]]]Exp[-t1];av[[1]]=fac(gbcdfC[r,bb[[1]]]-gbcdfC[l,bb[[1]]]);
 Do[z=u[[j]];fac=polyC[t1-z[[1]],aa[[j]]]Exp[z[[1]]-t1];av[[j]]+=fac(gbcdfC[r-z[[2]],bb[[j]]]-gbcdfC[l-z[[2]],bb[[j]]]);
  z=u[[j+ng-1]];fac=polyC[t1-z[[1]],bb[[j]]]Exp[z[[1]]-t1];av[[j]]+=fac(gbcdfC[r-z[[2]],aa[[j]]]-gbcdfC[l-z[[2]],aa[[j]]]),{j,2,ng}];{cost,av}];
sp[t1_?NumericQ,l_?NumericQ,r_?NumericQ,amp_,aa_List,bb_List,ord_Integer] := Module[{x,wt,h,cost=0.,av=ConstantArray[0.,ng]},
 If[r<=l,Return[{cost,av}]];{x,wt}=quad[ord,l,r];Do[h=profile[likeAmp[x[[i]],amp]];cost+=wt[[i]]Exp[tau-t1-x[[i]]]phi[h];av+=wt[[i]]h pv[t1,x[[i]],aa,bb],{i,Length[x]}];{cost,av}];
run[aa_List,bb_List,no_Integer,ord_Integer] := Module[{ww=srcScale w0,tx,tw,amp,roots,z,a,cost=0.,av=ConstantArray[0.,ng],start=AbsoluteTime[]},
 {tx,tw}=outerQuad[no];Do[amp=sourceAmp[tx[[i]],ww];roots=boundary[#,amp]&/@{1.,N[n],phiPrime};
  {z,a}=cp[tx[[i]],roots[[1]],roots[[2]],1.,aa,bb];cost+=tw[[i]]z;av+=tw[[i]]a;
  {z,a}=cp[tx[[i]],roots[[3]],hi,cap,aa,bb];cost+=tw[[i]]z;av+=tw[[i]]cap a;
  {z,a}=sp[tx[[i]],roots[[2]],roots[[3]],amp,aa,bb,ord];cost+=tw[[i]]z;av+=tw[[i]]a,{i,Length[tx]}];
 <|"scaled_cost"->cost,"actual_logJ"->Log[cost]-tau,"min_average"->Min[av/coeff],"min_group"->First@Ordering[av/coeff,1],"group24_average"->(av/coeff)[[24]],"seconds"->AbsoluteTime[]-start|>];
compact[x_Association]:=KeyTake[x,{"scaled_cost","actual_logJ","min_average","min_group","group24_average","seconds"}];

srcScale=1.0276924765700282;cg=c+.10;j=24;baseA=ConstantArray[cg,ng];baseB=ConstantArray[cg,ng];
one[d1_,d2_,no_:128,ord_:32]:=Module[{aa=baseA,bb=baseB},aa[[j]]+=d1;bb[[j]]+=d2;<|"d_c1"->d1,"d_c2"->d2,"result"->compact[run[aa,bb,no,ord]]|>];
base=one[0.,0.]; hs={.005,.01,.0125,.025,.05};probes=Flatten[Table[{one[h,0.],one[-h,0.],one[0.,h],one[0.,-h]},{h,hs}],1];
at[d1_,d2_]:=FirstCase[probes,x_/;x["d_c1"]==d1&&x["d_c2"]==d2];
grads=Table[<|"h"->h,"d_group24_dc1"->(at[h,0.]["result","group24_average"]-at[-h,0.]["result","group24_average"])/(2h),"d_group24_dc2"->(at[0.,h]["result","group24_average"]-at[0.,-h]["result","group24_average"])/(2h)|>,{h,hs}];
best=probes[[First@Ordering[probes[[All,"result","min_average"]],-1]]];
localValidation=If[best["result","min_average"]>base["result","min_average"]+10^-10,<|"N256_soft32"->one[best["d_c1"],best["d_c2"],256,32],"N256_soft64"->one[best["d_c1"],best["d_c2"],256,64]|>,<|"not_run"->True|>];

global[offset_]:=Module[{r=evaluateReceiverC[ConstantArray[cg+offset,ng],srcScale,128,32]},<|"offset"->offset,"receiver_c"->cg+offset,"result"->compact[Join[r,<|"group24_average"->r["averages_by_group"][[24]],"min_group"->First@Ordering[r["averages_by_group"],1]|>]]|>];
globalControls=global /@ {.005,.01,.0125,.025};globalRefine=global[.0075];
winner=Join[globalControls,{globalRefine}][[First@Ordering[Join[globalControls,{globalRefine}][[All,"result","min_average"]],-1]]];
ans=<|"status"->"r034 direct fixed-source receiver audit; floating numerical evidence only","seed_file"->"rounds/r030/receiver_seed.json",
 "input_dimension"->n,"block_size"->b,"source_scale"->srcScale,"r033_global_receiver_c"->cg,"group"->24,"mapping"->"u[j]:(c1,c2), u[j+32]:(c2,c1)",
 "baseline_N128_soft32"->base,"coordinate_probes_N128_soft32"->probes,"central_group24_gradients"->grads,"best_local_coordinate"->best,"local_256_validation"->localValidation,
 "uniform_global_controls_N128_soft32"->globalControls,"one_global_refinement_N128_soft32"->globalRefine,"best_uniform_global"->winner,
 "interpretation"->"The minimum can switch groups and is nonsmooth. Uniform global c is a necessary comparator; no FOSD or receiver-specific gain is claimed."|>;
Export[outFile34,ans,"RawJSON"];Print[<|"bestLocal"->best["result","min_average"],"bestGlobal"->winner["result","min_average"]|>];
