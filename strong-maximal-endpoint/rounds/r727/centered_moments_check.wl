m3=Integrate[(1-w)^3Exp[-w],{w,0,1}]+Integrate[(w-1)^3Exp[-w],{w,1,Infinity}];
m4=Integrate[(w-1)^4Exp[-w],{w,0,Infinity}];
Print[FullSimplify[m3]];
Print[m4];
