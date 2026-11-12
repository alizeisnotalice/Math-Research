q=1-e;
tests={FullSimplify[q-1/(1+d+y)>d/2,0<d<=1/10&&0<e<d/4&&0<y<d],FullSimplify[4(2/q-1)<=12,0<e<=1/2],FullSimplify[Integrate[p/u^2,{u,p,b}]==1-p/b,0<p<b<=1]};
Print[tests];If[!And@@tests,Exit[1]];
