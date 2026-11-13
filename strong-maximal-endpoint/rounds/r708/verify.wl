(* Exact decimal-rational witness inputs; independent 60/100 digit physical scores. *)
folder=DirectoryName[$InputFileName];data=Import[FileNameJoin[{folder,"steps.json"}],"RawJSON"];
parse[s_String]:=Module[{a=StringSplit[s,"."]},If[Length[a]==1,FromDigits[First[a]],FromDigits[a[[1]]<>a[[2]]]/10^StringLength[a[[2]]]]];
heights=parse/@data["heights"];tailheight=parse[data["tail_height"]];tstart=parse[data["t0"]];tend=data["T"];panels=data["panels"];n=data["n"];k=n-1;
physical[depth_,prec_]:=Module[{hs=N[heights,prec],last=N[tailheight,prec],t0=N[tstart,prec],dt=N[(tend-tstart)/panels,prec],H=N[depth,prec],aa,boundaries,values,cost},
 aa[x_?NumericQ]:=Fold[1+#2 #1/x&,1,Range[k]];
 boundaries=Table[t0+j dt,{j,0,panels}];
 values=(Exp[t0-#](1+#/H)^k aa[H+#]&)/@boundaries;
 cost=Total[MapThread[#1(1+Log[#1])^k #2&,{hs,Most[values]-Rest[values]}]]+last(1+Log[last])^k Last[values];
 <|"depth"->depth,"precision"->prec,"log_score"->ToString[N[t0+Log[aa[H]]-Log[cost],prec-10],InputForm]|>
];
values=Flatten[Table[physical[H,prec],{H,{4194304,1073741824}},{prec,{60,100}}]];
mu=1+ProductLog[1/E];ps=(mu-1)/mu;
mm[p_]:=1+p Log[(1-p)/p];
second=FullSimplify[D[mm[p],{p,2}]+1/(p(1-p)^2),Assumptions->0<p<1];
curvature=Factor[(1-2p)^2/2+2p(1-p)^2-1/2];
response[p_]:=1+Log[(1-p)/(1-q)]+u p;
stationary=FullSimplify[D[response[p],p]/.{p->q,u->1/(1-q)}];
shape=FullSimplify[D[response[p],{p,2}],Assumptions->0<p<1];
out=<|"version"->$Version,"physical_scores"->values,"m_second_identity"->ToString[second,InputForm],
 "response_stationary_identity"->ToString[stationary,InputForm],"response_second"->ToString[shape,InputForm],
 "sign_ratio_quadratic"->ToString[curvature,InputForm],"positive_curvature_magnitude"->ToString[N[2ps^2(1-ps),90],InputForm]|>;
Export[FileNameJoin[{folder,"wolfram.json"}],out,"RawJSON"];Print[out];Quit[];
