folder=DirectoryName[$InputFileName];
Clear[bb,tt,mm,nn,ll,zz,W];ff=Exp[-bb^2/(2mm tt)]/Sqrt[2Pi mm tt];psi=ff(nn/tt)^ll W[bb/Sqrt[mm tt]];
raw=FullSimplify[(D[psi,tt]-mm D[psi,{bb,2}]/2)/(ff(nn/tt)^ll/tt),Assumptions->mm>0&&tt>0&&nn>0&&Element[ll,Reals]];
identity=FullSimplify[(raw/.bb->Sqrt[mm tt]zz)-(-ll W[zz]+zz W'[zz]/2-W''[zz]/2),Assumptions->mm>0&&tt>0&&nn>0&&Element[{zz,ll},Reals]];
w=2-zz-Exp[zz];ou=FullSimplify[(D[w,{zz,2}]-zz D[w,zz])/2];
calc[job_]:=Module[{rec=job[[1]],prec=job[[2]],t,u,wall,mu,sig,c,p,K,L,n0,C2,M,b,f0,pot,v0,int,res,off},
 t=rec["t"];u=rec["u"];wall=TrueQ[rec["wall"]];mu=N[1+ProductLog[1/E],prec+25];sig=Sqrt[mu];
 c=N[Rationalize[rec["c"],0],prec+25];p=1/(1+Exp[c]);off=Round[100rec["fee_offset"]]/100;K=(mu-1)/8+off;
 L=N[Log[t],prec+25];n0=t/2;C2=If[wall,6,0];M=If[wall,1,0];b=sig Sqrt[t]u;f0=Exp[-u^2/2]/(sig Sqrt[2Pi t]);
 pot[at_,xx_?NumericQ]:=Module[{z=xx/(sig Sqrt[at]),f},f=Exp[-z^2/2]/(sig Sqrt[2Pi at]);Erfc[-z/Sqrt[2]]/2-(N[(9-mu)/24,prec+20]z^2+C2+K Log[at/n0])f-M L f (n0/at)^(1/32)(2-z-Exp[z])];
 v0=pot[t-1,b];
 int=MemoryConstrained[TimeConstrained[NIntegrate[Exp[-x]((1-p)(pot[t-1,b+mu-x]-v0)+p(pot[t-1,b+mu-x-c]-v0)),{x,0,Infinity},WorkingPrecision->prec,AccuracyGoal->prec-20,PrecisionGoal->prec-20,MaxRecursion->20,Method->{"GlobalAdaptive","SymbolicProcessing"->0}],50,$Failed],268435456,$Failed];
 res=If[int===$Failed,$Failed,t*(int+v0-pot[t,b])/f0];
 <|"t"->t,"u"->u,"wall"->wall,"fee_offset"->rec["fee_offset"],"precision"->prec,"c"->ToString[c,InputForm],"scaled_residual"->ToString[res,InputForm],"python_residual"->rec["scaled_residual"]|>
];
records=Import[FileNameJoin[{folder,"residual.json"}],"RawJSON"];
selected=Select[records,(MemberQ[{16384,1048576},#["t"]]&&MemberQ[{-2,0},#["u"]]&&Abs[#["fee_offset"]-.02]<10^-12)||(#["t"]==1048576&&#["u"]==0&&!TrueQ[#["wall"]]&&#["fee_offset"]<-.01)||(MemberQ[{4096,16384,65536,262144},#["t"]]&&#["u"]==-4&&TrueQ[#["wall"]]&&Abs[#["fee_offset"]-.02]<10^-12)&];
jobs=Flatten[Table[{rec,prec},{rec,selected},{prec,{60,100}}],1];
LaunchKernels[Min[12,Length[jobs]]];DistributeDefinitions[calc];values=ParallelMap[calc,jobs,Method->"CoarsestGrained"];CloseKernels[];
values=DeleteDuplicatesBy[values,{#["t"],#["u"],#["wall"],#["fee_offset"],#["precision"]}&];
out=<|"version"->$Version,"heat_identity"->ToString[identity,InputForm],"OU_expression"->ToString[ou,InputForm],
 "lambda"->"1/32","OU_compact_margin"->"7/96","OU_left_margin_coefficient"->"13/32","quadrature"->values|>;
Export[FileNameJoin[{folder,"wolfram.json"}],out,"RawJSON"];Print[<|"heat_identity"->identity,"samples"->Length[values],"failed"->Count[values,r_/;r["scaled_residual"]=="$Failed"]|>];Quit[];
