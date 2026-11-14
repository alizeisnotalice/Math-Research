(* Algebraic parameter derivative: q'(c)=-q(1-q). *)
op[x_] := D[x, c] - q (1-q) D[x, q];
mm = 1+c q;
vv = 1+c^2 q (1-q);
kk = 2+c^3 q (1-q) (1-2q);
Print[Factor[op[op[mm]] /. {c -> 1+r, q -> r/(1+r)}]];
Print[Factor[op[-kk/(6vv)] /. {c -> 1+r, q -> r/(1+r)}]];
