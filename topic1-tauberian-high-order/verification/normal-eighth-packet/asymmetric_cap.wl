ClearAll[eps,z,x,y,p,q,f,a,h,cut];
Print["Wolfram: ",$Version];
Print["Model: 0<eps<1/10, 0<=z<=1. p=1-eps z, q=1+eps z for plus; reversed for minus; f=1/eps."];
a=p/(f-p);h=p^2/f;cut=a(p-q);
w=2q-x-y;
plus=Integrate[Integrate[w,{y,0,p-x/a}],{x,0,h}]+Integrate[Integrate[w,{y,0,a(p-x)}],{x,h,p}];
minus=Integrate[Integrate[w,{y,0,q}],{x,0,cut}]+Integrate[Integrate[w,{y,0,p-x/a}],{x,cut,h}]+Integrate[Integrate[w,{y,0,a(p-x)}],{x,h,q}];
pl=FullSimplify[plus/.{p->1-eps z,q->1+eps z,f->1/eps},0<eps<1/10&&0<=z<=1];
mi=FullSimplify[minus/.{p->1+eps z,q->1-eps z,f->1/eps},0<eps<1/10&&0<=z<=1];
Print["plus series: ",FullSimplify[Normal[Series[pl,{eps,0,2}]]]];
Print["minus series: ",FullSimplify[Normal[Series[mi,{eps,0,2}]]]];
Print["sum series: ",FullSimplify[Normal[Series[pl+mi,{eps,0,2}]]]];
Print["sum quadratic z coefficient: ",Coefficient[Coefficient[Normal[Series[pl+mi,{eps,0,2}]],eps,2],z,1]//FullSimplify];
(* Check split order 0<=cut<=h<=q<p for minus. *)
pt=1+eps z;qt=1-eps z;ft=1/eps;ht=pt^2/ft;ct=pt/(ft-pt)(pt-qt);
Print["minus branch guard violations: ",Reduce[0<eps<1/10&&0<=z<=1&&(ct<0||ct>ht||ht>qt||qt>pt),{eps,z},Reals]];
