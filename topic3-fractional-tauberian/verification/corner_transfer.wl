ClearAll[x,y,a,k];b[z_]:=Max[0,-a,z-a/2];dd=a+b[x]+b[y];m=Min[x,y];
gap=Resolve[ForAll[{x,y,a},Implies[x>=0&&y>=0,dd-x-y>=Max[0,Abs[a]-2m]/2]],Reals];
u[z_]:=Min[k,Max[0,k+z+a/2]]+Min[k,Max[0,k-z+a/2]];
h=a*dd+b[x]*b[y]-u[x]*u[y]+2k^2;
g=a^2/4+a*(x+y)/2+x*y+2k^2-(k+Max[0,k-x+a/2])*(k+Max[0,k-y+a/2]);
ident=Resolve[ForAll[{x,y,a,k},Implies[x>=0&&y>=0&&k>=0&&Abs[a]<=2m,h==g]],Reals];
Print[<|"gap"->gap,"plateau_identity"->ident,"version"->$Version|>];
