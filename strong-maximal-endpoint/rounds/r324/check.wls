a=3/25;rate=Log[2]+Log[21/19]+(21/10)a/(1-1/40);checks={rate<11/10,11/10<1+ProductLog[1/E]};Print[{checks,N[rate,20]}];If[!And@@checks,Exit[1]];
