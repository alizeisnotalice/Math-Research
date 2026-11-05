(* Exact shell-integral, parameter, limit, and power-certificate audit. *)
ClearAll[guard,a,s,n,gamma,cc];
guard[name_,test_]:=If[TrueQ[test],Print["PASS ",name],Print["FAIL ",name];Exit[1]];
records=Flatten[Table[
 Module[{nn=n,gg=gamma,c0=cc,beta,d,tau,e,P,p,tailPower,pcoef,m1,m2,A,r0,radialI,uShell,fShell,delta},
  beta=1-gg;d=nn gg c0;tau=nn beta+d;e=nn-tau;P=d/gg;p=nn gg/e;tailPower=d/e;pcoef=(nn-1)gg/d;
  m1=tau gg/e;m2=tau(nn-1)/((d/gg)*(nn-gg));A=tau/(nn-gg)pcoef^gg;
  r0=a^(gg/e);radialI=Integrate[s^(P-1),{s,r0,1},Assumptions->a>0];
  uShell=a^p+(nn-1)a radialI; fShell=a^m1+tau(nn-1)/(nn-gg)a^beta radialI;
  delta=If[d/tau<1/3,d/tau,1/4];
  {TrueQ[FullSimplify[0<d&&0<e&&nn beta<tau<nn&&And@@Thread[{P,p,tailPower,pcoef,m1,m2,A}>0]]],
   TrueQ[FullSimplify[uShell==a^p+pcoef a(1-a^tailPower),a>0]],
   TrueQ[FullSimplify[fShell==a^m1+m2 a^beta(1-a^tailPower),a>0]],
   TrueQ[FullSimplify[Limit[fShell/uShell^beta,a->0,Direction->"FromAbove"]==A]],TrueQ[0<delta<1/3]}],
 {n,{2,3}},{gamma,{1/4,1/2,3/4,9/10}},{cc,{1/100,1/20,1/10,1/4,1/2,3/4,19/20}}],2];
guard["all_parameter_positivity",And@@(#[[1]]&/@records)];guard["all_shell_volume_identities",And@@(#[[2]]&/@records)];guard["all_shell_mass_identities",And@@(#[[3]]&/@records)];guard["all_leading_limits",And@@(#[[4]]&/@records)];guard["all_delta_policies_feasible",And@@(#[[5]]&/@records)];
asChecks=Flatten[Table[Module[{nn=n,gg=gamma,beta=1-gamma,d,tau,e,p,tailPower,pcoef,m1,m2,A,u,f},
 d=(nn-1)gg;tau=nn-gg;e=nn-tau;p=nn gg/e;tailPower=d/e;pcoef=(nn-1)gg/d;m1=tau gg/e;m2=tau(nn-1)/((d/gg)*(nn-gg));A=tau/(nn-gg)pcoef^gg;u=a^p+pcoef a(1-a^tailPower);f=a^m1+m2 a^beta(1-a^tailPower);TrueQ[FullSimplify[u==a&&f==a^beta&&A==1,a>0]]],{n,{2,3}},{gamma,{1/4,1/2,3/4,9/10}}],1];guard["AS_radial_zero_exact",And@@asChecks];
acChecks=Flatten[Table[Module[{nn=n,gg=gamma,beta=1-gamma,dv,tau,e,p,tailPower,pcoef,m1,m2,u,f,aa=1/16},
 tau=nn beta+dv;e=nn-tau;p=nn gg/e;tailPower=dv/e;pcoef=(nn-1)gg/dv;m1=tau gg/e;m2=tau(nn-1)/((dv/gg)*(nn-gg));u=aa^p+pcoef aa(1-aa^tailPower);f=aa^m1+m2 aa^beta(1-aa^tailPower);TrueQ[FullSimplify[Limit[u,dv->0,Direction->"FromAbove"]==aa(1+(nn-1)/nn Log[1/aa])]]&&TrueQ[FullSimplify[Limit[f,dv->0,Direction->"FromAbove"]==aa^beta(1+(nn-1)beta/(nn-gg)Log[1/aa])]]],{n,{2,3}},{gamma,{1/4,1/2,3/4,9/10}}],1];guard["AS2_d_to_zero_AC_limit",And@@acChecks];
B10=(2/3)^10 5^5/11;guard["special_delta",1/10/(11/10)==1/11];guard["special_B_power",B10==3200000/649539];guard["special_B10_gt_(7_over_6)^10",B10>(7/6)^10];Print["PASS round31 symbolic audit; parameter triples=",Length[records]];
