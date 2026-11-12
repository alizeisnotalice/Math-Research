nlist={9,4,8,3};r=3; lhs=Total[nlist^2]/r;rhs=Total[nlist]^2/(r Length[nlist]);tests={lhs>=rhs,FullSimplify[(1+t)-(1-a)==a+t]};Print[{lhs,rhs,tests}];If[!TrueQ[And@@tests],Exit[1]];
