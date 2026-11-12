s=1-1/60;r=1-33/1000;k=h(1-(r+s)/2);
checks={s>r,FullSimplify[s k+(1-s)(k-h)==h(s-r)/2],
FullSimplify[r k+(1-r)(k-h)==-h(s-r)/2],
Resolve[ForAll[{l,c,h},Implies[l>0&&c>0&&h>2l(c+1)/(s-r),h(s-r)/(2l)-c>1]],Reals]};
Print[checks];If[!And@@checks,Exit[1]];
