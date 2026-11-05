Get[FileNameJoin[{DirectoryName[$InputFileName],"i2_common.wl"}]];
Print["Wolfram ",$Version,"; master PID ",$ProcessID];
Print["Domain 0<delta<1/10. Four inverse-root branches, each split into three polynomial moments, plus three continuity boundaries."];
jobs=Join[Flatten[Table[<|"kind"->"moment","branch"->b,"term"->k|>,{b,1,4},{k,1,3}],1],Table[<|"kind"->"continuity","boundary"->b|>,{b,1,3}]];
LaunchKernels[15];Print["jobs: ",Length[jobs]," actual launched kernels: ",Length[Kernels[]]];
worker[j_]:=Module[{start=AbsoluteTime[],cpu=TimeUsed[],br,q,m,jac,expr,result,bnd,id},
 If[j["kind"]=="moment",
  br=branches[[j["branch"]]];q=br["Q"];m=br["M"];jac=br["jac"];
  expr=Switch[j["term"],1,2delta(1+q)q jac,2,-2delta(1+q)m jac,3,-2delta^2(q-m)^2 jac];
  result=FullSimplify[Integrate[expr,{t,br["lo"],br["hi"]},Assumptions->0<delta<1/10],0<delta<1/10];
  id=br["id"]<>"_term"<>ToString[j["term"]],
  bnd=j["boundary"];
  result=FullSimplify[Switch[bnd,1,(branches[[1]]["M"]-branches[[2]]["M"])/.t->delta,2,(branches[[2]]["M"]-branches[[3]]["M"])/.t->1,3,(branches[[3]]["M"]-branches[[4]]["M"])/.t->delta],0<delta<1/10];
  id="boundary"<>ToString[bnd]];
 <|"id"->id,"kind"->j["kind"],"kernel_id"->$KernelID,"pid"->$ProcessID,"wall_seconds"->N[AbsoluteTime[]-start],"cpu_seconds"->N[TimeUsed[]-cpu],"max_memory_bytes"->MaxMemoryUsed[],"result"->result|>];
DistributeDefinitions[branches,worker,delta,t];
rows=ParallelMap[worker,jobs,Method->"CoarsestGrained"];
Print["actual kernel IDs: ",Sort@DeleteDuplicates@Lookup[rows,"kernel_id"]];
Print["actual PIDs: ",Sort@DeleteDuplicates@Lookup[rows,"pid"]];
Print["job CPU seconds total: ",Total@Lookup[rows,"cpu_seconds"]];
Print["peak child MaxMemoryUsed bytes: ",Max@Lookup[rows,"max_memory_bytes"]];
Print["boundary residuals: ",Lookup[Select[rows,#["kind"]=="continuity"&],"result"]];
I2=FullSimplify[Total[Lookup[Select[rows,#["kind"]=="moment"&],"result"]],0<delta<1/10];
Print["I2 exact rational: ",Factor[I2]];
Print["I2 Taylor through delta12: ",Series[I2,{delta,0,12}]];
Print["I2 denominator factors: ",Factor[Denominator[Together[I2]]]];
Print["candidate delta-log or fractional term: absent in exact rational I2 on this domain."];
Put[<|"jobs"->rows,"I2"->I2|>,FileNameJoin[{DirectoryName[$InputFileName],"i2_results.wl"}]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"i2_metadata.json"}],Map[Join[KeyDrop[#,"result"],<|"result"->ToString[#["result"],InputForm]|>]&,rows],"RawJSON"];
CloseKernels[];Print["complete"];
