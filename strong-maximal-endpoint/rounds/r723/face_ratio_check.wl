ph = (1-q+q Exp[I c t])/(1-I t);
mm = 1+c q;
ps = Exp[-I mm t] ph;
rr = Exp[I c t]/ph;
Print[Simplify[(D[rr,t]/.t->0)-I(c-mm)]];
Print[Factor[(D[rr,{t,2}]/.t->0)/.{c->1+r,q->r/(1+r)}]];
Print[Simplify[(D[Log[ps],{t,2}]/.t->0)+(1+c^2 q(1-q))]];
