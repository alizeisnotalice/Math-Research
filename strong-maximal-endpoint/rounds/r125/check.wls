p=1/2;v=1+8Log[2];f[a_]=2Log[p/a]-(1+v)(p-a)+4Log[1/p]-2(1-p);s=1/(1+4Log[2]);
checks={FullSimplify[f'[s]==0],FullSimplify[f''[a]>0,0<a<1/2],FullSimplify[f[1/4]==4Log[2]-3/2]};
Print[{checks,N[{s,f[s]},20]}];If[!And@@checks,Exit[1]];
