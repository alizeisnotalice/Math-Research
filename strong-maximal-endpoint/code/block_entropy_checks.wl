base=DirectoryName[DirectoryName[$InputFileName]];
Clear[t,c];mu=1+ProductLog[Exp[-1]];
i2[prec_]:=Module[{cc=N[mu,prec+30],pp,rr,B,pieces},rr=Exp[-cc];pp=rr/(1+rr);B=1+Exp[cc];pieces={
 NIntegrate[Exp[-t]t(2pp)^2/(1+rr)^2,{t,0,cc},WorkingPrecision->prec,AccuracyGoal->prec-25,PrecisionGoal->prec-25],
 NIntegrate[Exp[-t](3t-2cc)(2pp(1-B/(3t-2cc)))^2/(1+rr)^2,{t,cc,2cc},WorkingPrecision->prec,AccuracyGoal->prec-25,PrecisionGoal->prec-25],
 NIntegrate[Exp[-t](4t-4cc)(2pp(1-2B/(4t-4cc)))^2/(1+rr)^2,{t,2cc,Infinity},WorkingPrecision->prec,AccuracyGoal->prec-25,PrecisionGoal->prec-25]};Total[pieces]];
p[c_]:=1/(1+Exp[c]);m[c_]:=1+c p[c];v[c_]:=1+c^2 p[c](1-p[c]);
mp=Together[D[m[c],c]];vp=Together[D[v[c],c]];
mpResidual=FullSimplify[(mp/.c->x)/.Exp[x]->1/(x-1)];vpResidual=FullSimplify[((vp/.c->x)/.Exp[x]->1/(x-1))-(x-1)];
mppResidual=FullSimplify[((Together[D[m[c],{c,2}]]/.c->x)/.Exp[x]->1/(x-1))+(x-1)/x];
normalMoment=Integrate[(z^2-1)^2 Exp[-z^2/2]/Sqrt[2Pi],{z,-Infinity,Infinity}];
k3=2+c^3 r(1-r)(1-2r);k4=6+c^4 r(1-r)(1-6r(1-r));
k3Residual=FullSimplify[(k3/.r->(c-1)/c)-c(3-c)];
k4Residual=FullSimplify[(k4/.r->(c-1)/c)-c(3-c)(4-c)];
gramResidual=FullSimplify[(c(3-c)(4-c)-(c(3-c))^2/c)-c(3-c)];
result=<|"wolfram"->$Version,"I2_60digits"->ToString[N[i2[60],35],InputForm],"I2_90digits"->ToString[N[i2[90],65],InputForm],"mean_derivative_critical_residual"->ToString[mpResidual,InputForm],"second_mean_derivative_critical_residual"->ToString[mppResidual,InputForm],"variance_derivative_critical_residual"->ToString[vpResidual,InputForm],"normal_centered_square_moment"->normalMoment,"critical_cumulant_3_residual"->k3Residual,"critical_cumulant_4_residual"->k4Residual,"critical_Gram_residual"->gramResidual,"stationary_fisher_limit_70digits"->ToString[N[p[mu]^2/2,70],InputForm]|>;

symbolicResult=result;
Clear[t];mu=1+ProductLog[Exp[-1]];
blockI[bb_Integer,prec_Integer]:=Module[{cc=N[mu,prec+30],r,p,B,poly,shift,parts},r=Exp[-cc];p=r/(1+r);B=1+Exp[cc];parts=Table[
 poly=Sum[Binomial[bb,k](t-k cc)^(bb-1),{k,0,j}];shift=If[j==0,0,Sum[Binomial[bb-1,k](t-(k+1)cc)^(bb-2),{k,0,j-1}]];
 NIntegrate[Exp[-t]poly/(1+r)^bb/Gamma[bb]*(bb p(1-B(bb-1)shift/poly))^2,{t,j cc,If[j<bb,(j+1)cc,Infinity]},WorkingPrecision->prec,AccuracyGoal->prec-15,PrecisionGoal->prec-15],{j,0,bb}];Total[parts]];
rows=Table[With[{val=blockI[bb,50]},Print[{bb,N[val,18]}];<|"b"->bb,"working_precision"->50,"fisher_35digits"->ToString[N[val,35],InputForm]|>],{bb,{2,4,8,16}}];
check=blockI[4,80];result=<|"wolfram"->$Version,"status"->"piecewise high-precision quadrature; not interval certified","rows"->rows,"I4_80precision_60digits"->ToString[N[check,60],InputForm]|>;

fisherResult=result;
seed=Import[FileNameJoin[{base,"rounds","r030","receiver_seed.json"}],"RawJSON"];
run=Import[FileNameJoin[{base,"rounds","r031","checks.json"}],"RawJSON"]["receiver_repair"];
rat[x_]:=Rationalize[SetPrecision[x,Infinity],0];c=rat[seed["c"]];b=seed["n"]/2;cap=seed["height_cap"];tau=rat[seed["tau_ref"]];lo=rat[run["finite_domain"][[1]]];hi=rat[run["finite_domain"][[2]]];vv=rat/@seed["v_nonnegative"];dd=rat/@seed["receiver_d"];
weights=Table[Binomial[b,k]Exp[-c k]/(1+Exp[-c])^b,{k,0,b}];
tail[z_]:=Total[Table[weights[[k+1]] GammaRegularized[b,Max[z-c k,0]],{k,0,b}]];
tails=Table[With[{u1=dd[[j]]+vv[[j]],u2=dd[[j]]-vv[[j]]},With[{t1=tail[hi-u1],t2=tail[hi-u2]},cap(t1+t2-t1 t2)]],{j,Length[vv]}];
nuTail=2Exp[-lo-hi]-Exp[-2hi];phi=cap(1+Log[cap])^(seed["n"]-1);ct=phi nuTail;
result=<|"wolfram"->$Version,"status"->"80-digit reevaluation of stored binary inputs; no directed interval certificate","normalization_residual"->ToString[(Total[N[weights,130]]-1),InputForm],"min_average_tail_80digits"->ToString[N[Min[tails],80],InputForm],"max_average_tail_80digits"->ToString[N[Max[tails],80],InputForm],"scaled_cost_tail_80digits"->ToString[N[Exp[tau]ct,80],InputForm],"actual_cost_tail_80digits"->ToString[N[ct,80],InputForm]|>;

tailResult=result;

result=<|"symbolic"->symbolicResult,"high_precision_fisher"->fisherResult,"receiver_tail"->tailResult|>; Export[FileNameJoin[{base,"rounds","r031","symbolic_checks.json"}],result,"RawJSON"]; Print[result];
