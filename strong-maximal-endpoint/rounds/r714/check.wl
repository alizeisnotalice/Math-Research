folder=DirectoryName[$InputFileName];
verify[M_,a_,b_]:=Module[{boxes,ends,cells,mass,inc,vol,x,cons,primal,lam,z,t,subsets,hall,dual,resp,allsub,subsetvals},
 boxes={{0,M},{-a,b},{M-b,M+a}};ends=Union[Flatten[boxes]];cells=Partition[ends,2,1];mass=Subtract@@Reverse[#]&/@cells;
 inc=Table[Boole[box[[1]]<=cell[[1]]&&cell[[2]]<=box[[2]]],{box,boxes},{cell,cells}];vol=Subtract@@Reverse[#]&/@boxes;x=Array[xx,Length[mass]];
 cons=Join[Thread[(inc.DiagonalMatrix[mass]).x>=vol/2],Thread[x>=0],Thread[x<=1]];
 primal=First[Minimize[{mass.x,And@@cons},x]];
 lam=Array[ll,3];z=Array[zz,Length[mass]];t=vol.lam;subsets=Rest[Subsets[Range[3]]];
 (* Explicit atom support union: Boole threads elementwise. *)
 hall=Table[Total[(vol lam)[[ss]]]<=Total[Table[If[Total[inc[[ss,j]]]>0,mass[[j]],0],{j,Length[mass]}]]/Total[mass]*t,{ss,subsets}];
 dual=t/2-mass.z;resp=First[Maximize[{dual,And@@Join[Thread[lam>=0],Thread[z>=0],Thread[z>=Transpose[inc].lam-1],hall]},Join[lam,z]]];
 subsetvals=Table[Module[{rows=inc[[ss]],v=vol[[ss]],value,uv},value=First[Minimize[{mass.x,And@@Join[Thread[(rows.DiagonalMatrix[mass]).x>=v/2],Thread[x>=0],Thread[x<=1]]},x]];uv=Total[Table[If[Total[rows[[All,j]]]>0,mass[[j]],0],{j,Length[mass]}]];<|"rows"->ss,"source"->ToString[value,InputForm],"gain"->ToString[uv/value,InputForm]|>],{ss,subsets}];
 <|"M"->M,"a"->a,"b"->b,"source"->ToString[primal,InputForm],"responsible_testing"->ToString[resp,InputForm],"true_formula_pass"->(primal==M/2),"responsible_formula_pass"->(resp==(a+b)(M+2a)/(2(b+2a))),"relative_loss"->ToString[1-resp/primal,InputForm],"subcatalogues"->subsetvals|>];
params=Join[{{10,1,3}},Select[Tuples[{{120},{1,5,11,15,20},{30,36,42,48,54}}],#[[2]]+#[[3]]<60&]];
LaunchKernels[12];DistributeDefinitions[verify];res=ParallelMap[Function[p,TimeConstrained[MemoryConstrained[verify@@p,268435456,$Failed],120,$Aborted]],params];CloseKernels[];
opt=Maximize[{x(1/2-2x)/(1/4+2x),0<x<1/4},x];out=<|"version"->$Version,"results"->res,"family_relative_loss_maximum"->ToString[opt,InputForm],"loss_y_derivative"->ToString[Factor[D[x(1-2x-2y)/(y+2x),y]],InputForm],"maximum_identity"->FullSimplify[First[opt]-(2-Sqrt[3])/4]|>;
Export[FileNameJoin[{folder,"wolfram.json"}],out,"RawJSON"];Print[<|"cases"->Length[res],"formula_passes"->Count[res,z_Association/;TrueQ[z["true_formula_pass"]&&z["responsible_formula_pass"]]]|>];Quit[];
