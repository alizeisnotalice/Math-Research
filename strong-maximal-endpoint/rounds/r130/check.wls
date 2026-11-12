den=1+d+y;b=1/den;
as=0<d<=1/10&&0<y<d&&1-3d/5<=q<=1-2d/5&&1<x<=4/3;
tests={FullSimplify[(1+y)/den<q&&b<q,as],FullSimplify[2(1-q)-(x-1)(q-b)>=2d/15,as],FullSimplify[2(4/3-(2/q-1))>=58/141,as],Limit[(1/(2(1-3d/5))-1/(2(1-2d/5)))/d,d->0]==1/10};Print[tests];If[!And@@tests,Exit[1]];
