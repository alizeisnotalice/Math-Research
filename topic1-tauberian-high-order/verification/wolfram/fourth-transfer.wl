(* Formal transfer conditional on proving C_n=L_n+O_n(a v^4), n>=3.
   These computations do not establish that analytic theorem. *)
Clear[v,n,m,u,am,expr,poly,rows,failures,record];
Print["VERSION: ",$Version];LaunchKernels[15];Print["KERNELS: ",Length[Kernels[]]];
rows=ParallelTable[Module[{p,got,want},
 p=Normal[Series[1/(1-v^n/(1+v)^n),{v,0,n+3}]];
 got=Table[Coefficient[p,v,n+j],{j,0,3}];
 want={1,-n,n(n+1)/2,-Binomial[n+2,3]+KroneckerDelta[n,3]};
 {n,got,got===want,$KernelID}],{n,3,17},Method->"CoarsestGrained"];
Scan[Print["linear_coefficients: ",ToString[#,InputForm]]&,rows];
failures=Count[rows,{_,_,False,_}];
SetAttributes[record,HoldRest];record[name_,test_]:=Module[{r=test},Print[name,": ",ToString[r,InputForm]];If[r=!=True,failures++]];
am=(m-1)/(2m);
expr=1+8(u^3+am u^6)-48 u^4(1+am u^3)^(4/3)+192 u^5(1+am u^3)^(5/3)-576(u^3+am u^6)^2;
poly=Normal[Series[expr,{u,0,6}]];
record["n3_uniform_m_fourth_coefficient",FullSimplify[Coefficient[poly,u,6]==-572-4/m,m>=1]];
record["n_ge4_linear_parameter_error_absorbed",FullSimplify[2>=1+4/n,n>=4]];
record["n_ge3_fractional_parameter_error_absorbed",FullSimplify[2+1/n>=1+4/n,n>=3]];
record["uniform_parameter_cubic_error_absorbed",FullSimplify[3>=1+4/n,n>=3]];
Print["Failures: ",failures];CloseKernels[];Exit[If[failures==0,0,1]];
