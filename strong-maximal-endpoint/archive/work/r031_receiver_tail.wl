seed=Import["/Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal/strong-maximal-endpoint/rounds/r030/receiver_seed.json","RawJSON"];
run=Import["/Users/zhengzhihao/Desktop/T/work/sm_receiver_continuous31.json","RawJSON"];
rat[x_]:=Rationalize[SetPrecision[x,Infinity],0];c=rat[seed["c"]];b=seed["n"]/2;cap=seed["height_cap"];tau=rat[seed["tau_ref"]];lo=rat[run["finite_domain"][[1]]];hi=rat[run["finite_domain"][[2]]];vv=rat/@seed["v_nonnegative"];dd=rat/@seed["receiver_d"];
weights=Table[Binomial[b,k]Exp[-c k]/(1+Exp[-c])^b,{k,0,b}];
tail[z_]:=Total[Table[weights[[k+1]] GammaRegularized[b,Max[z-c k,0]],{k,0,b}]];
tails=Table[With[{u1=dd[[j]]+vv[[j]],u2=dd[[j]]-vv[[j]]},With[{t1=tail[hi-u1],t2=tail[hi-u2]},cap(t1+t2-t1 t2)]],{j,Length[vv]}];
nuTail=2Exp[-lo-hi]-Exp[-2hi];phi=cap(1+Log[cap])^(seed["n"]-1);ct=phi nuTail;
result=<|"wolfram"->$Version,"status"->"80-digit reevaluation of stored binary inputs; no directed interval certificate","normalization_residual"->ToString[(Total[N[weights,130]]-1),InputForm],"min_average_tail_80digits"->ToString[N[Min[tails],80],InputForm],"max_average_tail_80digits"->ToString[N[Max[tails],80],InputForm],"scaled_cost_tail_80digits"->ToString[N[Exp[tau]ct,80],InputForm],"actual_cost_tail_80digits"->ToString[N[ct,80],InputForm]|>;
Export["/Users/zhengzhihao/Desktop/T/work/r031_receiver_tail.json",result,"RawJSON"];Print[result];
