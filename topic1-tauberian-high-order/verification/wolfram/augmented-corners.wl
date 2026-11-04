(* Exact scalar certificates for the proposed all-critical-coordinate argument.
   Geometric capture and all-input transfer require separate proofs. *)
Clear[qq,xx,vv,aa,bb,dd,y,check,results,failures];
Print["VERSION: ",$Version];
LaunchKernels[15];Print["KERNELS: ",Length[Kernels[]]];
results=ParallelTable[Module[{poly,certificate,derivative,rootvalue},
 poly=xx^qq-qq xx+qq-1;
 certificate=(xx-1)^2 Sum[(qq-1-j)xx^j,{j,0,qq-2}];
 {qq,Expand[poly-certificate]===0,
  And@@Table[qq-1-j>0,{j,0,qq-2}],
  Simplify[(qq+qq/(qq-1)-4)==(qq-2)^2/(qq-1)],
  $KernelID}],{qq,2,16},Method->"CoarsestGrained"];
Scan[Print["q_certificate: ",ToString[#,InputForm]]&,results];
check[name_,expr_]:=Module[{r=expr},Print[name,": ",ToString[r,InputForm]];If[r=!=True,failures++]];
failures=Count[results,{_,False,_,_,_}|{_,_,False,_,_}|{_,_,_,False,_}];
check["all_q_exponent_at_least_four",
 Reduce[qq>=2 && qq+qq/(qq-1)<4,qq,Reals]===False];
check["bounded_positive_clipping",
 Resolve[ForAll[{xx,aa,bb},Implies[0<=xx<=aa && 0<bb<=aa,
    Min[xx,bb]>=bb xx/aa]],Reals]];
Print["Failures: ",failures];CloseKernels[];
Exit[If[failures==0,0,1]];
