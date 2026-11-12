h[z_]=2/z-1;p=1/2;a=1/4;
d1=Integrate[h[z]-h[p],{z,a,p}];d2=2a/p Integrate[h[z],{z,p,1}];
checks=FullSimplify[{d1==2Log[2]-1,d2==2Log[2]-1/2,d1+d2==4Log[2]-3/2}];
Print[{checks,N[d1+d2,30]}];If[!And@@checks,Exit[1]];
