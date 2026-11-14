base=DirectoryName[$InputFileName];
f[b_,s_]:=Exp[b]((1-s)-Exp[-s b])/s;
series=Normal[Series[f[m+a s+c s^2,s],{s,0,2}]];
(* exp(m)(m-1)=1, m=mu. *)
sol1=FullSimplify[Coefficient[series,s,1]/.a->m/2];
sol2=FullSimplify[Coefficient[series,s,2]/.{a->m/2,c->m(9-m)/24}];
k=Normal[Series[t(m+(m/2)(t/(1+t))+(m(9-m)/24)(t/(1+t))^2),{t,0,3}]];
mu=1+ProductLog[Exp[-1]];r=Exp[-mu];
psi[s_]:=-Log[1-s]+Log[(1+r Exp[mu s])/(1+r)];
res=<|"b1_residual"->sol1,"b2_residual"->sol2,"K_series"->ToString[Expand[k],InputForm],"psi_variance_at_zero"->ToString[FullSimplify[D[psi[s],{s,2}]/.s->0/.Exp[-mu]->mu-1],InputForm],"uniform_third_derivative_lower"->ToString[N[2-mu^3/(6Sqrt[3]),70],InputForm],"gamma_c0"->ToString[N[Sqrt[2Pi mu](1/2-CDF[NormalDistribution[],2]/3-2CDF[NormalDistribution[],-1]/3),70],InputForm]|>;
res["three_point_mean_residue"]=Expand[(h+a)/2+(1/2-pp)(h-a)+pp b-(h-pp(h-a-b))];
res["three_point_mgf_residue"]=FullSimplify[(Exp[t(h+a)]/2+(1/2-pp)Exp[t(h-a)]+pp Exp[t b])Exp[-t mm]-(Exp[-dd t]Cosh[a t]-pp Exp[-(dd+a)t]+pp Exp[t(b-mm)])/.h->mm-dd];
res["three_point_variance_residue"]=Expand[(h+a-(h-pp(h-a-b)))^2/2+(1/2-pp)(h-a-(h-pp(h-a-b)))^2+pp(b-(h-pp(h-a-b)))^2-(a^2+pp((h-b)^2-a^2)-pp^2(h-a-b)^2)];
res["truncated_exponential_half_moment_identity"]=FullSimplify[2(1-z)/(1-z^2)-2/(1+z),Assumptions->0<=z<1];
res["truncated_exponential_half_moment_violation"]=ToString[Reduce[0<=z<=1&&2/(1+z)>2,z,Reals],InputForm];
Export[FileNameJoin[{base,"lp_expansion_wolfram.json"}],res,"RawJSON"];Print[res];Quit[];
