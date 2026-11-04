(* Formal coefficients conditional on the candidate a(v^6+v^(n+1)) upper.
   This script does not establish that all-input geometric theorem. *)
Clear[n,p,v,m,u,am,rows,record,failures];Print["VERSION: ",$Version];
LaunchKernels[15];Print["KERNELS: ",Length[Kernels[]]];
rows=ParallelTable[Module[{p=Min[6,n+1],s,got,want,last},
 s=Normal[Series[1/(1-v^n/(1+v)^n),{v,0,n+p-1}]];
 got=Table[Coefficient[s,v,n+j],{j,0,p-1}];
 want=Table[(-1)^j Binomial[n+j-1,j]+KroneckerDelta[j,n],{j,0,p-1}];
 {n,p,got,got===want,Last[got]2^(n+p-1),$KernelID}],{n,2,16},Method->"CoarsestGrained"];
Scan[Print["linear_coefficients: ",ToString[#,InputForm]]&,rows];failures=Count[rows,{_,_,_,False,_,_}];
SetAttributes[record,HoldRest];record[label_,test_]:=Module[{r=test},Print[label,": ",ToString[r,InputForm]];If[r=!=True,failures++]];
am=(m-1)/(2m);
Do[Module[{s=u^n(1+am u^n),root=u(1+am u^n)^(1/n),coefficient,expected},
 coefficient=Coefficient[Normal[Series[1/(1-2^n s/(1+2root)^n),{u,0,2n}]],u,2n];
 expected=2^(2n)((-1)^n Binomial[2n-1,n]+1)+2^(n-1)(1-1/m);
 Print["quadratic_parameter_coefficient_n",n,": ",ToString[FullSimplify[coefficient,m>=1],InputForm]];
 record["uniform_m_coefficient_n"<>ToString[n],FullSimplify[coefficient==expected,m>=1]]],{n,2,5}];
record["n_ge6_parameter_quadratic_absorbed",FullSimplify[2>=1+6/n,n>=6]];
record["fractional_parameter_error_absorbed_low_dimensions",FullSimplify[2+1/n>=1+(n+1)/n,n>=2]];
record["fractional_parameter_error_absorbed_high_dimensions",FullSimplify[2+1/n>=1+6/n,n>=6]];
Print["FAILURES: ",failures];CloseKernels[];Exit[If[failures==0,0,1]];
