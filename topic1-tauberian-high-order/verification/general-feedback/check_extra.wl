ClearAll[x,y,z,q,a];
f[x_,y_]:=Max[x,0]Max[y,0](Max[x,0]+Max[y,0]);
Print["convexity clipped violation:",Reduce[x>=0&&y>=0&&-Max[x,y]<=z<=Max[x,y]&&(f[x+z,y+z]+f[x-z,y-z])/2<x y(x+y),{x,y,z},Reals]];
Print["powers factored:",FullSimplify[q^a>=1,Assumptions->q>=1&&a>=2]];
