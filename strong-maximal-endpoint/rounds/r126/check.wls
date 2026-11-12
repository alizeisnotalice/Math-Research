p=1/2;h[z_]=2/z-1;v=1+8Log[2];q[a_]=(2Log[p/a]-(p-a)+4Log[2]-1)/(p-a);
t0=t/.FindRoot[t-Log[t]==1/2+2Log[2],{t,3},WorkingPrecision->60];
a0=p/t0;s=1/(1+4Log[2]);
checks={Abs[q[a0]-h[a0]]<10^-40,a0<s,Abs[q'[a0]]<10^-40,q[a0]>v};
Print[{checks,N[{a0,s,q[a0],v},15]}];If[!And@@checks,Exit[1]];
