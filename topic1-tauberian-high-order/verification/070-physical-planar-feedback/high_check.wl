ClearAll[del,b,t,aa,x];
P[z_]:=z^2/(2del)+2aa z^3/3;
qhi=1-t/del-aa t^2+t;mhi=2P[b]-P[t]+t^2/2;
Print["Wolfram ",$Version,"; PID ",$ProcessID];
Print["join with Phi(b)=1: ",Factor[((qhi-b)/.t->b)/.(aa->(del-b)/(del b^2))]];
Print["high q-M under Phi(b)=1: ",Factor[(qhi-mhi)/.(aa->(del-b)/(del b^2))]];
Print["high q-M by t=b x: ",Factor[(qhi-mhi)/.{aa->(del-b)/(del b^2),t->b x}]];
Print["nonnegative violation: ",Reduce[0<b<del<1/8&&0<=(del-b)/(del b^2)<1/4&&0<=x<=1&&((qhi-mhi)/.{aa->(del-b)/(del b^2),t->b x})<0,{del,b,x},Reals]];
