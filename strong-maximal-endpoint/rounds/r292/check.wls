Clear[p,a,c,r,l,v,h];
checks={
Resolve[ForAll[{p,a,v,l,h},Implies[0<p<1&&a>1&&v>0&&0<h<l<6/p&&h>6(a+1)v/p^2,(p+a v)(1-h/l)<p-v]],Reals],
Resolve[ForAll[{c,r,v,p},Implies[c>1&&0<r<1/(2c)&&v>0,((p+v)+r(p-c v))/(1+r)>p]],Reals]};
Print[checks];If[checks=!={True,True},Exit[1]];
