Clear[t];mu=1+ProductLog[Exp[-1]];
blockI[bb_Integer,prec_Integer]:=Module[{cc=N[mu,prec+30],r,p,B,poly,shift,parts},r=Exp[-cc];p=r/(1+r);B=1+Exp[cc];parts=Table[
 poly=Sum[Binomial[bb,k](t-k cc)^(bb-1),{k,0,j}];shift=If[j==0,0,Sum[Binomial[bb-1,k](t-(k+1)cc)^(bb-2),{k,0,j-1}]];
 NIntegrate[Exp[-t]poly/(1+r)^bb/Gamma[bb]*(bb p(1-B(bb-1)shift/poly))^2,{t,j cc,If[j<bb,(j+1)cc,Infinity]},WorkingPrecision->prec,AccuracyGoal->prec-15,PrecisionGoal->prec-15],{j,0,bb}];Total[parts]];
rows=Table[With[{val=blockI[bb,50]},Print[{bb,N[val,18]}];<|"b"->bb,"working_precision"->50,"fisher_35digits"->ToString[N[val,35],InputForm]|>],{bb,{2,4,8,16}}];
check=blockI[4,80];result=<|"wolfram"->$Version,"status"->"piecewise high-precision quadrature; not interval certified","rows"->rows,"I4_80precision_60digits"->ToString[N[check,60],InputForm]|>;
Export["/Users/zhengzhihao/Desktop/T/work/r031_wolfram_fisher.json",result,"RawJSON"];Print[result];
