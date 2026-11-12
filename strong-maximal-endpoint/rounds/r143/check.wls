q=(p-e)/(1-e); tests={FullSimplify[q+(1-q)e==p,0<e<p<1],9/20+(1-9/20)/10==101/200,((1/2-1/10)/(1-1/10))==4/9,101/200-1/2==1/200};
Print[tests];If[!TrueQ[And@@tests],Exit[1]];
