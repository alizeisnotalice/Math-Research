ClearAll[q];
Print["Mathematica: ",$Version];
beta=q(2q-1)/(q-1);
vals={beta+2-19/2,q+q/(q-1)+5-19/2,(q+4)+2-9};
Print["factors: ",Factor /@ vals];
Print["violations q>=3: ",Reduce[q>=3&&#<0,q,Reals]& /@ vals];
Print["boundary q=3: ",vals/.q->3];
