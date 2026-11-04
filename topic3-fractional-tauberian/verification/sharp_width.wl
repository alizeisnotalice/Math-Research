ClearAll[t,n,b,p,l,s,al];
F=(l+s)/2+s^(1-n)(1-al(l^p+s^p));h=l+l^(1-n)-2al*l^(1-b);
Print[<|"F_l"->FullSimplify[D[F,l]],"F_s_at_unit"->FullSimplify[D[F,s]/.{s->1,l->1,al->1,p->n-b}],"h_prime_at_unit"->FullSimplify[D[h,l]/.{l->1,al->1}],"h_second_central"->FullSimplify[D[h,{l,2}]/.{l->1,al->1,b->n/2}],"two_dimensional_width_series"->ToString[Series[1/(1-t)-(1-t),{t,0,8}],InputForm],"geometric_remainder"->FullSimplify[(1/(1-t)-(1-t))-(2t+Sum[t^j,{j,2,6}])-t^7/(1-t)],"n3_half_series"->ToString[Series[(1-t)^(-4/3)-(1-t)^(2/3),{t,0,5}],InputForm],"version"->$Version|>];

Print[FullSimplify[(2^(-b/n))^(n/b),Assumptions->n>0&&b>0]];
