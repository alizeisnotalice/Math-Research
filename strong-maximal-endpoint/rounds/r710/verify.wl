folder=DirectoryName[$InputFileName];
parse[s_String]:=Module[{q=StringSplit[s,"/"]},If[Length[q]==1,FromDigits[q[[1]]],FromDigits[q[[1]]]/FromDigits[q[[2]]]]];
lpcheck[rec_]:=Module[{boxes=rec["boxes"],n=rec["n"],axes,tuples,lo,hi,mat,keep,mass,vol,x,res,value,interval},
 axes=Table[Union[Flatten[boxes[[All,j]]],If[TrueQ[rec["common_point"]],{0},{}]],{j,n}];
 tuples=Tuples[Range[Length[#]-1]&/@axes];lo=Table[Table[axes[[j,q[[j]]]],{j,n}],{q,tuples}];hi=Table[Table[axes[[j,q[[j]]+1]],{j,n}],{q,tuples}];
 (* Explicit indexing avoids any pairing ambiguity in the Cartesian atomization. *)
 mat=Table[Boole[And@@Table[boxes[[i,j,1]]<=lo[[a,j]]&&hi[[a,j]]<=boxes[[i,j,2]],{j,n}]],{i,Length[boxes]},{a,Length[tuples]}];
 keep=Select[Range[Length[tuples]],Total[mat[[All,#]]]>0&];mat=mat[[All,keep]];lo=lo[[keep]];hi=hi[[keep]];
 mass=Times@@@ (hi-lo);vol=Times@@@(Subtract@@@Reverse[#,2]&/@boxes);
 x=Array[xx,Length[mass]];
 res=MemoryConstrained[TimeConstrained[Minimize[{mass.x,And@@Join[Thread[(mat.DiagonalMatrix[mass]).x>=vol/2],Thread[x>=0],Thread[x<=1]]},x],60,$Failed],536870912,$Failed];
 value=If[res===$Failed,$Failed,First[res]];interval=parse/@rec["plain"]["source_interval"];
 <|"n"->n,"seed"->rec["seed"],"common_point"->rec["common_point"],"atoms"->Length[mass],"exact_minimum"->ToString[value,InputForm],
 "inside_python_certificate"->If[NumberQ[value],interval[[1]]<=value<=interval[[2]],False]|>
];
corecheck[rec_]:=Module[{h=rec["h"],den=FromDigits[rec["a_denominator"]],num=FromDigits[rec["a_numerator"]],a,ratio,upper,exact},
 a=num/den;exact=(1+h^2a/(1+a))/(h+Log[2,h]);ratio=N[exact,100];
 <|"h"->h,"row_volume_inequality"->((den+num)^h<=2den^h<(den+num+1)^h),
 "exact_ratio_matches"->(exact==parse[rec["exact_gap_over_log_n"]]),"ratio100"->ToString[ratio,InputForm]|>
];
records=Import[FileNameJoin[{folder,"catalogues.json"}],"RawJSON"];
chosen=Select[records,#["n"]<=3&&#["K"]==2&];
cores=Import[FileNameJoin[{folder,"core_family.json"}],"RawJSON"];
LaunchKernels[Min[12,Length[chosen]]];DistributeDefinitions[parse,lpcheck];lp=ParallelMap[lpcheck,chosen,Method->"CoarsestGrained"];CloseKernels[];
core=corecheck/@cores;
r=Exp[-c];u0=(1-r)/(1+r);qlow=-Log[1-(1+r)u];qhigh=Log[2/((1+r)(1-u))];
derivatives={D[qlow,c],D[qhigh,c]};
lowbound=FullSimplify[0<=u rr/(1-(1+rr)u)<=1,Assumptions->0<rr<1&&0<u<=(1-rr)/(1+rr)];
out=<|"version"->$Version,"LP_checks"->lp,"core_checks"->core,
 "quantile_derivatives"->ToString[derivatives,InputForm],"lower_branch_Lipschitz"->ToString[lowbound,InputForm],
 "entropy_discretization_loss"->ToString[N[Log[2]-Log[1-Exp[-1]],80],InputForm],"core_gap_limit"->ToString[N[Log[2],100],InputForm]|>;
Export[FileNameJoin[{folder,"wolfram.json"}],out,"RawJSON"];Print[<|"LP_cases"->Length[lp],"LP_passed"->Count[lp,r_/;TrueQ[r["inside_python_certificate"]]],"core_cases"->Length[core]|>];Quit[];
