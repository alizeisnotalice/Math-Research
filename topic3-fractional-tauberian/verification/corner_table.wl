ClearAll[x,k,u];f[u_]:=u^2+u+x(1-x)+2k^2-(k+Max[0,k-x+u])(k+Max[0,k-1+x+u]);
v={f[-x],f[x],f[Clip[x-k,{-x,x}]],f[Clip[1-x-k,{-x,x}]]};
rows={{0<=k<=1/2&&0<=x<=(1-k)/2,2x},{0<=k<=1/2&&(1-k)/2<=x<=1/2,2x+2k(1-k-2x)},{1/2<=k<=2/3&&0<=x<=(1-k)/2,2x},{1/2<=k<=2/3&&(1-k)/2<=x<=1/2,2((1-k)^2+(2k-1)x)},{2/3<=k<=1&&0<=x<=(1-k)/2,2x},{2/3<=k<=1&&(1-k)/2<=x<=1-k,2((1-k)^2+(2k-1)x)},{2/3<=k<=1&&1-k<=x<=k/2,2k*x},{2/3<=k<=1&&k/2<=x<=1/2,k(2k-1)+2(1-k)x},{k>=1&&0<=x<=1/2,2k(1-k)+2(2k-1)x}};
checks=Map[TimeConstrained[Resolve[ForAll[{x,k},Implies[#[[1]],Max@@v==#[[2]]]],Reals],30,"timeout"]&,rows];Print[checks];
