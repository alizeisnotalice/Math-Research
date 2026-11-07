(* Exact algebra and high-precision corroboration for terminal profile stability. *)
Clear[x,nn,theta,cc,mm];
cher=FullSimplify[-(x/(8nn))x+4nn(x/(8nn))^2==-x^2/(16nn),Assumptions->nn>0];
hder=FullSimplify[D[Log[1+Log[1+x]]-x/(1+x),x] == (x-Log[1+x])/((1+x)^2(1+Log[1+x])),Assumptions->x>0];
mstat=FullSimplify[(1+Exp[mm](1-mm))/(1+Exp[mm])^2 /. Exp[mm]->1/(mm-1)]==0;
mu=1+ProductLog[1/E];
mgf[tt_,c_]:=(1+Exp[-(1-tt)c])/((1+Exp[-c])(1-tt));
rows=Flatten[Table[Module[{res=mu tt+4tt^2-Log[mgf[tt,c]]},<|"theta"->ToString[tt,InputForm],"c"->ToString[c,InputForm],"log_mgf_margin_60"->ToString[N[res,60],InputForm],"pass"->TrueQ[N[res,80]>=0]|>],{tt,{1/10^6,1/100,1/10,1/2}},{c,{0,mu/2,mu,2,5,20}}],1];
If[!TrueQ[cher&&hder&&mstat]||!And@@Lookup[rows,"pass"],Exit[2]];
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/stability61.json"];
Export[out,<|"wolfram"->$Version,"chernoff_algebra"->cher,"threshold_derivative"->hder,"mean_stationarity"->mstat,"mgf_grid"->rows,"scope"->"Symbolic identities and high-precision MGF corroboration; uniform policy proof is analytic, no random sampling."|>,"RawJSON"];Print[{cher,hder,mstat,Length[rows]}];
