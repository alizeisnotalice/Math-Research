Clear[t,h,l,r,p,e,b,F];
checks={
 Resolve[ForAll[{t,h,l,r},Implies[0<=l<h<t<=r,(h-l)/(r-l)<=h/t]],Reals],
 Resolve[ForAll[{t,h,p,e,b,F},Implies[0<h<t && 0<e<p<1 && 0<b<=h/t && p<b F+(1-b)(p-e),F>p+e(t-h)/h]],Reals],
 Resolve[ForAll[{t,h,l,p,F},Implies[h>0 && 0<t<l<h/p && 0<p<1 && p<F(1-t/l),F-p>p^2 t/h]],Reals]};
Print[checks];If[checks=!={True,True,True},Exit[1]];
