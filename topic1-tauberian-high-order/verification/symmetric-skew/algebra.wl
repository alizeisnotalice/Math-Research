Clear[n,e,k,u,p,T,rows,flags];Print["VERSION: ",$Version];LaunchKernels[15];Print["KERNELS: ",Length[Kernels[]]];
rows=ParallelTable[Module[{mat,det,expected,bal},
 mat=(1-e)IdentityMatrix[n]+e ConstantArray[1,{n,n}];det=Factor[Det[mat]];
 expected=(1-e)^(n-1)(1+(n-1)e);
 bal=Sum[Binomial[n,q]Binomial[q,q/2],{q,2,n,2}];
 {n,Simplify[det==expected],Coefficient[Expand[1-det],e,2]==Binomial[n,2],bal<3^n,
  FullSimplify[Coefficient[Expand[1-det],e,1]==0],$KernelID}],{n,2,16},Method->"CoarsestGrained"];
Scan[Print["dimension: ",ToString[#,InputForm]]&,rows];
p=(k^2-k^3/(8n))/256;
flags={FullSimplify[Binomial[n,2] (p/.k->16n/3)==Binomial[n,2]n^2/27,n>=2],
 FullSimplify[D[Binomial[n,2]p,k]==Binomial[n,2]k(16n-3k)/(2048n),n>=2],
 FullSimplify[(n+9)-(n+8)==1],
 Reduce[n>=2 && e>0 && (n-1)e<1/12 && 2(1-(n-1)e)<=1,{n,e},Reals]===False};
Print["general_scalar_checks: ",ToString[flags,InputForm]];
failures=Count[rows,{_,False,___}|{_,_,False,___}|{_,_,_,False,___}|{_,_,_,_,False,_}]+Count[flags,Except[True]];
Print["FAILURES: ",failures];CloseKernels[];Exit[If[failures==0,0,1]];
