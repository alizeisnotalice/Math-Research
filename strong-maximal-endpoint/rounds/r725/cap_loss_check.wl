i0=Integrate[(w-a)Exp[-w],{w,a,Infinity},Assumptions->Element[a,Reals]&&a>0];
i1=Integrate[(w+c-a)Exp[-w],{w,a-c,Infinity},Assumptions->Element[{a,c},Reals]&&a>0&&a>c];
Print[i0];
Print[i1];
Print[FullSimplify[((1-q)i0+q i1-2(1-q)Exp[-a])/.q->1/(1+Exp[c]),Assumptions->Element[{a,c},Reals]]];
