(* Actual matched reference payoff and Poisson terminal term. No interval certificate. *)
mu44=N[1+ProductLog[1/E],180];g44=(mu44-1)/8;cQ44=(9-mu44)/24;cR44=(3+mu44)/12;
args44=$ScriptCommandLine;out44=If[Length[args44]>=2,args44[[2]],"/tmp/endpoint44.json"];wp44=100;
mean44[a_]:=1+a/(1+Exp[a]);var44[a_]:=1+a^2 Exp[a]/(1+Exp[a])^2;
mgf44[a_,e_]:=-Log[1-e]+Log[(1+Exp[(e-1)a])/(1+Exp[-a])];y44[a_,e_]:=1/(1-e)+a Exp[(e-1)a]/(1+Exp[(e-1)a]);
z44[n_]:=z44[n]=zz/.FindRoot[Erfc[zz/Sqrt[2]]/2==Sqrt[mu44/n]zz,{zz,Sqrt[Log[n]]},WorkingPrecision->(wp44+30),AccuracyGoal->(wp44+10),PrecisionGoal->(wp44+10)];
params44[n_,co_]:=Module[{ll=Log[n],dd=Sqrt[mu44 n]z44[n],r,a,cc},cc=If[Precision[co]===Infinity,co,Rationalize[co,0]];r=n mu44-dd+cc ll;a=aa/.FindRoot[2aa-mean44[aa]==r/n,{aa,mu44-dd/(2n)},WorkingPrecision->(wp44+20),AccuracyGoal->(wp44-5),PrecisionGoal->(wp44-5)];{r,a,mean44[a]-a}];
exact44[n_Integer,co_?NumericQ]:=Module[{r,a,b,ks,xx,bi,ww,f,rr,cdf,sf,aux1,aux2,j,p,x2},
 {r,a,b}=params44[n,co];ks=Range[0,n];xx=Clip[r-a ks,{0,Infinity}];bi=N[Binomial[n,#]&/@ks,wp44+10];ww=bi Exp[-a ks]/(1+Exp[-a])^n;p=1/(1+Exp[a]);
 f=Exp[-r]Dot[bi,xx^(n-1)]/(Gamma[n](1+Exp[-a])^n);rr=Exp[-r]Dot[bi,xx^n]/(Gamma[n+1](1+Exp[-a])^n);cdf=Dot[ww,GammaRegularized[n,0,xx]];sf=Dot[ww,(r-a ks)GammaRegularized[n,0,xx]-n GammaRegularized[n+1,0,xx]];
 aux1=Exp[-r]Dot[bi,xx^(n+1)]/(Gamma[n+2](1+Exp[-a])^n);x2=Clip[r-a-a Range[0,n-2],{0,Infinity}];aux2=Exp[-r+a]Dot[N[Binomial[n-2,#]&/@Range[0,n-2],wp44+10],x2^(n-1)]/(Gamma[n](1+Exp[-a])^(n-2));j=p/8 (aux1+a^2 p(1-p)aux2)/f (EulerGamma+PolyGamma[0,N[n+1,wp44+20]]-1);
 <|"r"->r,"a"->a,"b"->b,"f"->f,"v_minus_one"->2b-cdf+rr+sf/n,"J_over_f"->j,"A_over_f"->(cdf-rr-sf/n)/f,"method"->"exact finite Gamma-binomial mixture"|>];
fourier44[n_Integer,co_?NumericQ,step_:1/8,cut_:24]:=Module[{r,a,b,lambda,e,y,vv,ss,ps,cn,fc,hA,aux,p,logpsi,center},
 {r,a,b}=params44[n,co];lambda=ll/.FindRoot[y44[a,-ll]==r/n,{ll,2b/var44[a]},WorkingPrecision->(wp44+20),AccuracyGoal->(wp44-5),PrecisionGoal->(wp44-5)];e=-lambda;y=r/n;vv=N[Range[-cut/step,cut/step]step,wp44];ss=e+I vv/Sqrt[n];p=1/(1+Exp[a]);
 logpsi=n(Log[(1-e)/(1-ss)]+Log[(1+Exp[(ss-1)a])/(1+Exp[(e-1)a])]-I y vv/Sqrt[n]);ps=Exp[logpsi];center=Re[Total[ps]]step/(2Pi Sqrt[n]);fc=Exp[n mgf44[a,e]-e r]center;
 hA=Re[Total[ps(-1/(ss(1-ss))-1/(n ss^2))]]/Re[Total[ps]];aux=Re[Total[ps(1/(1-ss)^2+a^2 p(1-p)Exp[a ss]/(1-p+p Exp[a ss])^2)]]/Re[Total[ps]];
 <|"r"->r,"a"->a,"b"->b,"lambda"->lambda,"f"->fc,"v_minus_one"->2b-fc hA,"J_over_f"->p/8 aux(EulerGamma+PolyGamma[0,N[n+1,wp44+20]]-1),"A_over_f"->hA,"method"->"negative tilted Fourier trapezoid","Fourier_step"->step,"Fourier_cutoff"->cut|>];
point44[n_,co_?NumericQ]:=If[n<=1024,exact44[n,co],fourier44[n,co,1/8,If[n>=2^512,32,24]]];
scaled44[n_,co_?NumericQ]:=Module[{d=point44[n,co]},d["v_minus_one"]/(d["f"]Log[n])];
json44[x_]:=N[x,18];
record44[n_,co_]:=Module[{d=point44[n,N[co,wp44+10]],l=Log[n]},<|"Fourier_cutoff"->If[n<=1024,0,If[n>=2^512,32,24]],"working_precision"->wp44,"N_exact"->ToString[n,InputForm],"C"->json44[co],"a"->json44[d["a"]],"r"->json44[d["r"]],"method"->d["method"],"density_N_over_logN"->json44[d["f"]n/l],"v_minus_one_over_f_logN"->json44[d["v_minus_one"]/(d["f"]l)],"J_over_f_logN"->json44[d["J_over_f"] /l],"reference_limit"->If[Abs[co-cQ44]<10^-50,0,json44[-(co-cQ44)]],"linearized_alpha_1_1_excess_over_f_logN"->json44[d["v_minus_one"]/(d["f"]l)+(11/10)d["J_over_f"]/l],"floor_Lambda1_over_f_logN"->json44[d["b"]^2/(d["f"]l)]|>];
If[!TrueQ[$endpoint44LibraryMode],
 rows44={};roots44={};dims44=Join[{2,4,8,16,32,64,128,256,1024},2^{12,18,26,40,64,96,128,256,512,1024}];
 Do[wp44=Max[120,2IntegerLength[n,10]+120];mu44=N[1+ProductLog[1/E],wp44+80];g44=(mu44-1)/8;cQ44=(9-mu44)/24;cR44=(3+mu44)/12;Do[row=record44[n,co];AppendTo[rows44,row],{co,{cQ44,cQ44+g44+1/10}}];Print[KeyTake[Last[rows44],{"N_exact","v_minus_one_over_f_logN","J_over_f_logN","linearized_alpha_1_1_excess_over_f_logN"}]];
 If[n>=16,cr=cc/.FindRoot[scaled44[n,cc]==0,{cc,N[cQ44,wp44],N[cQ44+1/10,wp44]},WorkingPrecision->70,AccuracyGoal->40,PrecisionGoal->40];AppendTo[roots44,<|"Fourier_cutoff"->If[n<=1024,0,If[n>=2^512,32,24]],"working_precision"->wp44,"N_exact"->ToString[n,InputForm],"reference_root_C"->json44[cr],"difference_from_cQ"->json44[cr-cQ44],"scaled_root_residual"->ToString[scaled44[n,cr],InputForm]|>]];
 Export[out44,<|"wolfram"->$Version,"working_precision_policy"->"Max[120,2 IntegerLength[N,10]+120]","status"->"actual matched payoff/root and exact Poisson diagnostics; not interval-certified","constants"-><|"mu"->json44[mu44],"g"->json44[g44],"cQ"->json44[cQ44],"cR"->json44[cR44],"target_alpha_1_1_margin_limit"->json44[-1/10+g44/10]|>,"rows"->rows44,"roots"->roots44|>,"RawJSON"],{n,dims44}];
];
