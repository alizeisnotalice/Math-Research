SetDirectory[DirectoryName[$InputFileName]];Print["VERSION ",$Version," PID ",$ProcessID];
(* Independent symbolic slab identities. *)
Fsym[s_]:=K(2 rr s+s^2);Jsym[s_]:=K(rr s^2+s^3/3);Hsym[s_]:=(2 ww-ee)s-s^2/2+Jsym[s]/ee;Csym[s_]:=Jsym[s]/ee-s^2/2;
Print["SYMBOLIC ",InputForm[{"Hprime"->FullSimplify[D[Hsym[z],z]==2ww-ee-z+Fsym[z]/ee],"Cprime"->FullSimplify[D[Csym[z],z]==Fsym[z]/ee-z],"rightSlab"->FullSimplify[Hsym[z]==2ww z-Integrate[ee+s-Fsym[s]/ee,{s,0,z}],ee>0]}]];
raw=Import["explicit.json","RawJSON"];u=2^-20;kap=16/3;lam=704/135;a=(2u)^2;w=u/(1+2u)+lam u^6;rho=1-2w;e=kap u^3;K=a;det=1-e^2;M=a rho^2;gapVol=(2w)^2-M;Nn=2^140;r=rho+w;
(* Rationalized exact positive root avoids subtracting near-equal radicals. *)
hh=N[(e^2/K)/(Sqrt[rho^2+e^2/K]+rho),250];FF[s_]:=K(2rho s+s^2);JJ[s_]:=K(rho s^2+s^3/3);
AA=2hh NIntegrate[e-FF[hh tt]/e,{tt,0,1},WorkingPrecision->235,AccuracyGoal->200,PrecisionGoal->200];
HH[s_]:=(2w-e)s-s^2/2+JJ[s]/e;CC[s_]:=JJ[s]/e-s^2/2;
z0=(gapVol-2AA)/(2w-e);zz=z0 (yy/.FindRoot[HH[z0 yy]==gapVol-2AA,{yy,1},WorkingPrecision->230,AccuracyGoal->190,PrecisionGoal->190]);cc=w-zz;
(* Direct nested integrals for the two bad corners, plus independent slab-cut integral. *)
bad1=2hh NIntegrate[1,{s,0,1},{b,0,e-FF[hh s]/e},WorkingPrecision->210,AccuracyGoal->175,PrecisionGoal->175];
bad2=2hh NIntegrate[1,{x,0,1},{y,0,e-FF[hh x]/e},WorkingPrecision->210,AccuracyGoal->175,PrecisionGoal->175];
cut=zz NIntegrate[FF[zz t]/e-zz t,{t,0,1},WorkingPrecision->210,AccuracyGoal->175,PrecisionGoal->175];
neighborDirect=2e^2-bad1-bad2-cut;neighborFormula=2e^2-2AA-CC[zz];defect=(M+neighborDirect)/(2r)^2;avg=1-defect;margin=u^2-defect;
Clear[decExact];decExact[s_String]:=Module[{p=StringSplit[s,RegularExpression["[eE]"]],b,x=0,sgn=1,dp},b=First[p];If[Length[p]>1,x=ToExpression[p[[2]]]];If[StringStartsQ[b,"-"],sgn=-1;b=StringDrop[b,1]];dp=If[StringContainsQ[b,"."],StringLength[Last[StringSplit[b,"."]]],0];sgn FromDigits[StringReplace[b,"."->""]]10^(x-dp)];
within[val_,key_]:=TrueQ[decExact[raw[key][[1]]]<=val<=decExact[raw[key][[2]]]];
checks={"AInside"->within[AA,"A"],"zInside"->within[zz,"z"],"cInside"->within[cc,"c"],"CInside"->within[CC[zz],"C"],"neighborInside"->within[neighborDirect,"neighbor_mass"],"meanInside"->within[avg,"average"],"marginInside"->within[margin,"strict_margin"],"cornerAgreement"->TrueQ[Abs[bad1-AA]<10^-180],"slabAgreement"->TrueQ[Abs[cut-CC[zz]]<10^-180],"neighborAgreement"->TrueQ[Abs[neighborDirect-neighborFormula]<10^-180],"strictMean"->TrueQ[margin>0]};
Print["CHECKS ",InputForm[checks]];Print["VALUES ",InputForm[{"hstar"->N[hh,35],"z"->N[zz,35],"c"->N[cc,35],"mean"->N[avg,35],"margin"->N[margin,35],"neighbor"->N[neighborDirect,35],"cornerError"->N[Abs[bad1-AA],5],"slabError"->N[Abs[cut-CC[zz]],5]}]];
