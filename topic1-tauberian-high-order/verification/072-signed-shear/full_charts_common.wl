ClearAll[e,w,f,s,y,z,x,term,chartIntegrand];
term[j_,d1_,d2_]:=Switch[j,1,e-d2,2,e-d1,3,(e-d1)(e-d2)/w];
chartIntegrand[chart_,j_]:=Switch[chart,
 1,Integrate[Integrate[Integrate[term[j,s,x],{x,s,e}],{s,0,e(e-y)/f}],{y,0,e}],
 3,Integrate[Integrate[Integrate[term[j,s+y,x],{x,s,e}],{s,0,e(e-y)/f}],{y,0,e}],
 2,Integrate[Integrate[Integrate[term[j,s+y,s],{s,0,(e-y)(e-z)/f}],{z,0,e}],{y,0,e}]];
