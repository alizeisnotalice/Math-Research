ClearAll[t,z,a,r,q,nu,eta2,M,g];a=1-t;r=t/2+t^2/2+31t^3/32;q=a+3t^2/16;nu=t^3/1024;eta2=t^3/16384;M=2q-a;
g[z_]:=(z-2r)^2-2a z+1;V[z_]:=4nu(z-2r)+4nu^2;
polys=Cancel/@{(g[a]-V[a]-4eta2)/t^3,(-g[q]-V[q])/t^3,a^2(-g[1/a]-V[1/a])/t^3,((a+2r)^2-4M(M-a)-1-16a nu-4eta2)/t^2};
low=Table[Coefficient[p,t,0]+Sum[Min[Coefficient[p,t,k],0]/100^k,{k,1,Exponent[p,t]}],{p,polys}];
Print[$Version];Print["patch polynomials:",polys];Print["whole interval lower:",low];If[!And@@Thread[low>0],Exit[1]];
rs=(a-Sqrt[2a^2-1])/2;
Print["normal depth bounds:",FullSimplify[t^3/32<=Sqrt[2](rs-r-nu) && Sqrt[2](rs-r+nu)<=t^3,0<t<=1/100]];
Print["patch separation:",FullSimplify[0<2r+2nu+2Sqrt[eta2]<a && a>=2r+4Sqrt[eta2] && a<q<(a+1/a)/2,0<t<=1/100]];
Print["transfer arithmetic:",{19(2^-24)^(2/3)<1/100,2Sqrt[2^-24]*(1/100)^(1/4)<1/2048,2^-24<2^-14}];

Print["moment coefficient:",FunctionExpand[Gamma[11/2]/2^24]];

ClearAll[cc,dd];ll=Sqrt[2Sqrt[2cc]dd+2dd^2];jj=((Sqrt[cc]+Sqrt[2]dd)ll-cc ArcSinh[ll/Sqrt[cc]])/2;
Print["sharp cap coefficient:",FullSimplify[Limit[jj/dd^(3/2),dd->0,Direction->"FromAbove"],cc>0]];
