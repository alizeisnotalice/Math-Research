r=4(6+y)/(4+y)-2;s=4(1+3/(4+y))-2;
gap=Integrate[2(s-r),{y,0,1}];
tests={FullSimplify[3/(4+y)>1/2&&s>r>3,0<y<1],FullSimplify[gap==8Log[5/4]],FullSimplify[(2(6+y)/(4+y))/(r+2)==1/2]};
Print[{tests,N[gap/14,18]}];If[!And@@tests,Exit[1]];
