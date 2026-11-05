Get[FileNameJoin[{DirectoryName[$InputFileName],"n2_common.wl"}]];
Print["Wolfram ",$Version," master PID ",$ProcessID];
Print["All 21 vertices; t∈[0,sqrt(2)-1], a=(1-t²)/(1+t²),b=2t/(1+t²), a>=b>=0."];
worker[{i_,j_}]:=Module[{f=area[i,j],limit=3/2-Sqrt[2]/4,bad,eq,cpu=TimeUsed[],start=AbsoluteTime[]},
 bad=Reduce[0<=t<=Sqrt[2]-1&&f>limit,t,Reals];
 eq=Reduce[0<=t<=Sqrt[2]-1&&f==limit,t,Reals];
 <|"i"->i,"j"->j,"formula"->ToString[Factor[f],InputForm],"bad"->ToString[bad,InputForm],"equality"->ToString[eq,InputForm],"kernel_id"->$KernelID,"pid"->$ProcessID,"cpu_seconds"->N[TimeUsed[]-cpu],"seconds"->N[AbsoluteTime[]-start],"max_memory_bytes"->MaxMemoryUsed[]|>];
jobs=Flatten[Table[{i,j},{i,1,6},{j,i,6}],1];
LaunchKernels[15];Print["jobs ",Length[jobs]," launched ",Length[Kernels[]]];
DistributeDefinitions[worker,area,vertex,la,lb,rat,t,aa,bb];
rows=ParallelMap[worker,jobs,Method->"FinestGrained"];
Print["actual IDs ",Sort@DeleteDuplicates@Lookup[rows,"kernel_id"]];Print["actual PIDs ",Sort@DeleteDuplicates@Lookup[rows,"pid"]];
Print["all violation sets False: ",And@@(#=="False"& /@ Lookup[rows,"bad"])];
Print["equality cases: ",Select[rows,#["equality"]!="False"&]];
Print["total CPU ",Total@Lookup[rows,"cpu_seconds"],"; max WL memory bytes ",Max@Lookup[rows,"max_memory_bytes"]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"n2_vertices.json"}],rows,"RawJSON"];
CloseKernels[];
