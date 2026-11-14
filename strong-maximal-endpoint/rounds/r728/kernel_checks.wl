ga=Exp[-x]x^(n-1)/Gamma[n];
Print[FullSimplify[D[ga,{x,2}]/ga-((x-(n-1))^2-(n-1))/x^2,Assumptions->x>0&&n>=3]];
op[f_]:=D[f,c]-q(1-q)D[f,q];
mm=1+c q;
vv=1+c^2 q(1-q);
Print[Factor[(2op[vv]^2/(-op[op[mm]]))/.{c->1+r,q->r/(1+r)}]];
