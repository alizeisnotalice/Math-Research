ClearAll[delta,q,d,s];
Print["Wolfram ",$Version,"; PID ",$ProcessID];
Print["Domain 0<delta<1/10, 0<=q,d<=1; H0(d)=delta(1-d)."];
t=(q-delta)/(1-delta);
mlo=FullSimplify[Integrate[q-d,{d,0,q}]+Integrate[q,{d,0,1-q/delta}]+Integrate[delta(1-d),{d,1-q/delta,1}],0<q<delta<1/10];
mhi=FullSimplify[Integrate[delta(1-d),{d,0,t}]+Integrate[q-d,{d,t,q}]+Integrate[delta(1-d),{d,0,1}],0<delta<q<1];
Print["M0 low q<=delta: ",Factor[mlo]];Print["M0 high q>=delta: ",Factor[mhi]];
Print["branch match: ",FullSimplify[(mlo-mhi)/.q->delta,0<delta<1/10]];
hlo=FullSimplify[delta(q-mlo),0<q<delta<1/10];hhi=FullSimplify[delta(q-mhi),0<delta<q<1];
Print["h1(q)=H1(1-q), low: ",Factor[hlo]];Print["h1(q)=H1(1-q), high: ",Factor[hhi]];
Print["h1 band violations low/high: ",{Reduce[0<delta<1/10&&0<=q<=delta&&(hlo<0||hlo>delta q),{delta,q},Reals],Reduce[0<delta<1/10&&delta<=q<=1&&(hhi<0||hhi>delta q),{delta,q},Reals]}];
(* Tonelli: in q-coordinate I(H)=2 integral [(1+q)h(q)-h(q)^2]. *)
I1=FullSimplify[2Integrate[(1+q)hlo-hlo^2,{q,0,delta}]+2Integrate[(1+q)hhi-hhi^2,{q,delta,1}],0<delta<1/10];
Print["I1 exact: ",Factor[I1]];Print["I1 Taylor through delta5: ",Series[I1,{delta,0,5}]];
I0=5delta/3-2delta^2/3;
Print["I0 minus I1: ",Factor[I0-I1]];
Print["I1 matches prior polynomial residual: ",Factor[I1-(I0-delta^2(155-155delta+86delta^2-22delta^3-4delta^4)/60)]];
Print["I1 rational in delta, no log/fractional powers by explicit form."];
