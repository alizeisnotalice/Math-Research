Clear[x,t,m,k,t0];h=k(x^2/(m t)+2Log[t/t0]);y=x-h;fy=Exp[-y^2/(2m t)]/Sqrt[2Pi m t];w=Erfc[-y/Sqrt[2m t]]/2;
exact=-k/t-k*x^2/(m*t^2)+2*k*x*h/(m*t^2)+2*k^2*x^2*(x-h)/(m^2*t^3);
res=FullSimplify[(D[w,t]-m D[w,{x,2}]/2)/fy-exact,Assumptions->{m>0,t>0,t0>0,Element[{x,k},Reals]}];
Print["nonlinear PDE residual = ",res];
p[v_]:=4+v^2+4v-10v^4+8v^5+v^6-12v^7/7;
Print["cap endpoint values = ", {p[0],p[1]}];
Print["cap derivatives matching = ", {p'[0]-4,p''[0]-2,p'''[0],p'[1],p''[1],p'''[1]}];
