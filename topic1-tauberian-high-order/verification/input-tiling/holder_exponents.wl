Print["VERSION ",$Version," PID ",$ProcessID];
zeta=q/(2q-1);delta=FullSimplify[q+6zeta-6];Print["factor ",InputForm[Factor[delta]]];
Print["factor identity ",InputForm[FullSimplify[delta==(2q-3)(q-2)/(2q-1),q>=2]]];
Print["q>=2 exponent negative ",InputForm[Reduce[q>=2 && delta<0,q,Reals]]];
Print["q>=3 exponent below 3/5 ",InputForm[Reduce[q>=3 && delta<3/5,q,Reals]]];
Print["q>=2 holder exponent below 1/3 ",InputForm[Reduce[q>=2 && 1-zeta<1/3,q,Reals]]];
Print["q3 equality ",InputForm[FullSimplify[delta/.q->3]]];
