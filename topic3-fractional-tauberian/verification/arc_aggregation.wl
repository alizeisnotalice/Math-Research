$Assumptions=Element[{t,zeta,d,e,h,k,n1,n3,s1,s3},Reals] && 0<t<1/100 && Abs[zeta]<=1/8;
ClearAll[t,zeta,d,e,z];a=1-t;c=2a^2-1;theta=Pi/4+zeta t;
rb=(a-Sqrt[c Tan[theta]])/2;sb=(a-Sqrt[c Cot[theta]])/2;
r=rb-d t^3 Cos[theta];s=sb-d t^3 Sin[theta];q=a+e t^2;M=2q-a;
g[z_]:=(z-2r)(z-2s)-2a z+1;
Print[$Version];Print["boundary coordinates:",FullSimplify[Normal@Series[#,{t,0,2}]]&/@{rb,sb}];
Print["leading margins:",FullSimplify[Limit[#/t^3,t->0,Direction->"FromAbove"]]&/@{g[a],-g[q],-g[1/a]}];
Print["mass leading:",FullSimplify[Limit[((a+2r)(a+2s)-4M(M-a)-1)/t^2,t->0,Direction->"FromAbove"]]];
Print["r,s boundary identity:",FullSimplify[(a+2rb)(a+2sb)-1-8rb sb,0<t<1/100 && -1/8<zeta<1/8]];
ClearAll[h,k,n1,n3,s1,s3];ux=Cos[theta];uy=Sin[theta];vx=-uy;vy=ux;
dx=2rb+(-2d+h(n1+n3))t^3 ux+k(s1+s3)t^(3/2)vx;
dy=2sb+(-2d+h(n1+n3))t^3 uy+k(s1+s3)t^(3/2)vy;
gp[z_]:=(z-dx)(z-dy)-2a z+1;
Print["pair perturbed leading:",FullSimplify[Limit[#/t^3,t->0,Direction->"FromAbove"]]&/@{gp[a],-gp[q],-gp[1/a]}];
Print["chosen parameter lower:",N[{2Sqrt[2](1/64-1/4096)-2/4096^2,1/4-2Sqrt[2](1/64+1/4096),2(1-1/64)-2Sqrt[2](1/64+1/4096),1-2/64},20]];

Print["relaxation homothety deficit:",Series[(2t^2-t^4)(1-a^2+c Log[c/a^2]),{t,0,6}]];
Print["quartic moment coefficient:",4!];

lowerBounds={2Sqrt[2](1/64-1/4096)-2/4096^2,1/4-2Sqrt[2](1/64+1/4096),2(1-1/64)-2Sqrt[2](1/64+1/4096),1-2/64};
If[!TrueQ[FullSimplify[And@@Thread[lowerBounds>0]]],Exit[1]];Print["exact leading lower bounds:",lowerBounds];
