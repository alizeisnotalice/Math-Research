args77=Rest[$ScriptCommandLine];src77=Import[args77[[1]],"RawJSON"];out77=args77[[2]];suffix77=If[StringContainsQ[args77[[1]],"_refine"],"_refine",""];If[!DirectoryQ[out77],CreateDirectory[out77,CreateIntermediateDirectories->True]];
ClearAll[rat77,log77,phi77,env77,rows77,case77];
rat77[s_String]:=If[StringMatchQ[s,RegularExpression["-?[0-9]+(/[0-9]+)?"]],ToExpression[s,InputForm],Throw["bad rational"]];
log77[t_]:=Module[{u=(t-1)/(t+1),lo,hi,den=10^12},lo=Sum[2u^(2k+1)/(2k+1),{k,0,19}];hi=lo+2u^41/(41(1-u^2));{Floor[den lo]/den,Ceiling[den hi]/den}];
phi77[t_]:=If[t<=1,{t,t},t(1+log77[t])];
env77[mesh_]:=Module[{nodes=Table[1+k/mesh,{k,0,mesh}],logs,low,up,a,b,va,vb,s},logs=log77/@nodes;low=Prepend[MapThread[{2+#2[[1]],-#1}&,{nodes,logs}],{1,0}];up={{1,0}};Do[a=nodes[[k]];b=nodes[[k+1]];va=a(1+logs[[k,2]]);vb=b(1+logs[[k+1,2]]);s=(vb-va)/(b-a);AppendTo[up,{s,va-s a}],{k,1,mesh}];{low,up}];
rows77[anchor_,cell_,ref_]:=Module[{step=2/ref,ans={},area,lx,rx,ly,ry,ox,oy},Do[area=Abs[(x-anchor[[1]])(y-anchor[[2]])];If[area!=0,{lx,rx}=Sort[{x,anchor[[1]]}];{ly,ry}=Sort[{y,anchor[[2]]}];ox=Table[Max[0,Min[rx,k+step]-Max[lx,k]],{k,0,8-step,step}];oy=Table[Max[0,Min[ry,k+step]-Max[ly,k]],{k,0,8-step,step}];AppendTo[ans,{Flatten[Outer[Times,ox,oy]],area}]],{x,cell[[1;;2]]},{y,cell[[3;;4]]}];ans];
case77[z_]:=Catch[Module[{ref=z["refinement"],eps=rat77[z["epsilon"]],sig=rat77[z["sigma"]],cells=z["receiver_cells"],f=rat77/@z["source"],w=rat77[z["source_cell_volume"]],mesh=z["mesh"],env,low,up,gs=z["groups"],used={},certs={},r,g,lam,theta,d,lb,ub,pc,pm,local,lower=0,upper=0,jlo,jhi,edges,expected,sh=z["shape"],af,lp,mat,rhs,cvec,sol,msgs,obj,independent=False,lpcount=0,t0=AbsoluteTime[]},
 env=env77[mesh];low=env[[1]];up=env[[2]];
 If[!And@@Thread[0<=f]||!And@@Thread[f<=2]||w!=1/ref^2||sig!=1+Min[eps/20,1/1000],Throw["source bounds"]];
 expected=Flatten[Table[If[sh=="notch",If[Quotient[i,ref]==0&&Quotient[j,ref]==0,1,1+eps],If[sh=="checker",If[EvenQ[Quotient[i,ref]+Quotient[j,ref]],2,1/2],If[MemberQ[{1,2},Quotient[i,ref]]||MemberQ[{1,2},Quotient[j,ref]],2,0]]],{i,0,4ref-1},{j,0,4ref-1}]];
 If[f!=expected,Throw["source definition"]];
 Do[
  used=Join[used,gg["receiver_cells"]];
  r=DeleteDuplicates[Flatten[rows77[gg["anchor"],cells[[#+1]],ref]&/@gg["receiver_cells"],1]];
  g=rat77/@gg["primal"];lam=ConstantArray[0,Length[r]];theta=ConstantArray[0,{Length[f],Length[low]}];
  Do[lam[[v[[1]]+1]]=rat77[v[[2]]],{v,gg["success_dual"]}];Do[theta[[v[[1]]+1,v[[2]]+1]]=rat77[v[[3]]],{v,gg["epigraph_dual"]}];
  If[Length[g]!=Length[f]||!And@@Thread[0<=g]||!And@@Thread[g<=f]||Min[r[[All,1]].g-sig r[[All,2]]]<0,Throw["primal constraints"]];
  If[Min[lam]<0||Min[Flatten[theta]]<0||!And@@Thread[(Total/@theta)<=w],Throw["dual constraints"]];
  d=theta.low[[All,1]]-Transpose[r[[All,1]]].lam;
  lb=sig(lam.r[[All,2]])+Total[Flatten[theta ConstantArray[low[[All,2]],Length[f]]]]+Total[MapThread[Min[0,#1 #2]&,{f,d}]];
  ub=w Total[(Max[up[[All,1]] #+up[[All,2]]]&)/@g];
  If[lb!=rat77[gg["lower"]]||ub!=rat77[gg["upper"]]||lb>ub,Throw["objective certificate"]];
  pc=Table[N[w Total[(If[#<=1,#,#(1+Log[#])]&)/@g],wp],{wp,{60,120}}];
  If[!And@@((TrueQ[N[lb,50]<=#<=N[ub,50]])&/@pc),Throw["numeric interval"]];
  (* An independent installed Wolfram LP solve for nontrivial small selected models. *)
  obj=Null;msgs={};
  If[ref==1&&z["assignment"]=="greedy"&&(sh=="checker"||sh=="cross"||(sh=="notch"&&eps==1/2)),
   mat=Join[(Join[#[[1]],ConstantArray[0,Length[f]]]&)/@r,
     Flatten[Table[Join[-low[[k,1]] UnitVector[Length[f],j],UnitVector[Length[f],j]],{j,Length[f]},{k,Length[low]}],1],
     (Join[-#,ConstantArray[0,Length[f]]]&)/@IdentityMatrix[Length[f]]];
   rhs=Join[sig r[[All,2]],Flatten[ConstantArray[low[[All,2]],Length[f]]],-f];cvec=Join[ConstantArray[0,Length[f]],ConstantArray[w,Length[f]]];
   If[Min[mat.Join[g,(Max[up[[All,1]] #+up[[All,2]]]&)/@g]-rhs]<0,Throw["matrix feasibility bug"]];
   Block[{$MessageList={}},sol=LinearProgramming[cvec,mat,rhs,Method->"Simplex"];msgs=ToString[#,InputForm]&/@$MessageList];
   If[!VectorQ[sol,NumericQ],Throw["Wolfram LP failed"]];If[!VectorQ[sol,(IntegerQ[#]||Head[#]===Rational)&]||Min[mat.sol-rhs]<0,Throw["LP not exact feasible"]];obj=cvec.sol;
   If[Abs[obj-gg["LP_objective"]]>10^-7||obj<N[lb,40]-10^-7||obj>N[ub,40]+10^-7,Throw["LP objective mismatch"]];lpcount++
  ];
  lower+=lb;upper+=ub;AppendTo[certs,<|"exact"->True,"rows"->Length[r],"primal_60"->ToString[pc[[1]],InputForm],"primal_120"->ToString[pc[[2]],InputForm],"Wolfram_LP_objective"->If[obj===Null,Null,ToString[N[obj,30],InputForm]],"Wolfram_LP_exact"->(obj=!=Null),"Wolfram_LP_messages"->msgs|>],
 {gg,gs}];
 If[Sort[used]!=z["covered_cells"]||Length[DeleteDuplicates[used]]!=Length[used],Throw["receiver partition"]];
 jlo=w Total[(First@phi77[#]&)/@f];jhi=w Total[(Last@phi77[#]&)/@f];
 If[jlo!=rat77[z["source_modular_lower"]]||jhi!=rat77[z["source_modular_upper"]]||lower!=rat77[z["allocation_lower"]]||upper!=rat77[z["allocation_upper"]],Throw["summed budget"]];
 If[Total[((cells[[#+1,2]]-cells[[#+1,1]])(cells[[#+1,4]]-cells[[#+1,3]])/4&)/@used]!=rat77[z["witness_area"]],Throw["receiver area"]];
 <|"shape"->sh,"epsilon"->z["epsilon"],"refinement"->ref,"assignment"->z["assignment"],"exact_pass"->True,"groups"->certs,"independent_Wolfram_LPs"->lpcount,"seconds"->N[AbsoluteTime[]-t0]|>
]];
sym77=TrueQ[FullSimplify[D[t(1+Log[t]),{t,2}]==1/t,Assumptions->t>1]];
LaunchKernels[12];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[rat77,log77,phi77,env77,rows77,case77];
t77=AbsoluteTime[];checks77=ParallelMap[TimeConstrained[case77[#],90,"TIMEOUT"]&,src77["rows"],Method->"FinestGrained"];CloseKernels[];
pass77=sym77&&AllTrue[checks77,AssociationQ[#]&&TrueQ[#["exact_pass"]]&];
Export[FileNameJoin[{out77,"verification77"<>suffix77<>".json"}],<|"version"->$Version,"workers"->12,"seconds"->N[AbsoluteTime[]-t77],"symbolic_phi_convexity"->sym77,"pass"->pass77,"rows"->checks77|>,"RawJSON"];
Print[If[pass77,"PASS ","FAIL "],Length[checks77]," independently reconstructed cases"];
If[!pass77,Print[Cases[checks77,_String]];Exit[3]];
