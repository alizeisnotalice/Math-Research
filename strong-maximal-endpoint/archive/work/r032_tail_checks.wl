Clear[t,s,n,mu];
inv2[u_]:=Exp[u]/ProductLog[Exp[u+1]];
invAsym=Limit[Exp[-t]t inv2[t-s],t->Infinity,Assumptions->Element[s,Reals]];
(* Fenchel coefficient with representative critical scale and subexponential epsilon. *)
coefLog=Log[n]+(n Log[2]-LogGamma[n]-(mu n-Sqrt[n Log[n]])+Log[n]^2)/(n-1);
coefLimit=Limit[coefLog,n->Infinity,Assumptions->mu>0];
mu0=1+ProductLog[Exp[-1]];
result=<|"wolfram"->$Version,"inverse_envelope_Gamma2_integrand_limit"->ToString[invAsym,InputForm],"Fenchel_log_coefficient_limit"->ToString[coefLimit,InputForm],"Fenchel_coefficient_limit_at_mu_50digits"->ToString[N[2Exp[1-mu0],50],InputForm],"fixed_N"->2,"s_for_numeric_check"->3,"Gamma2_inverse_integrand_samples"->Table[<|"t"->tt,"value_40digits"->ToString[N[Exp[-tt]tt inv2[tt-3],40],InputForm]|>,{tt,{20,100,1000}}]|>;
Export["/Users/zhengzhihao/Desktop/T/work/r032_tail_checks.json",result,"RawJSON"];Print[result];
