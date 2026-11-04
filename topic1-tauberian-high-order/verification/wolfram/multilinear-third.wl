Clear[m,eps,q,u,n,am,sigma,expr,s];
Print["Wolfram version: ",$Version];
failures=0;SetAttributes[record,HoldRest];
record[label_,expr_]:=Module[{r=expr},Print[label," = ",InputForm[r]];If[r=!=True,failures++]];
sigma=m(1-(1-eps)^(1/m));
record["Uniform quadratic coefficient",
 FullSimplify[(D[sigma,{eps,2}]/.eps->0)/2==(m-1)/(2m),m>=1]];
record["Third derivative identity",
 FullSimplify[D[sigma,{eps,3}]==(1-1/m)(2-1/m)(1-eps)^(1/m-3),
   m>=1 && 0<eps<1]];
record["Uniform derivative coefficient bound",
 Resolve[ForAll[q,Implies[0<=q<=1,0<=(1-q)(2-q)<=2]],Reals]];
record["Logarithmic form of the positive-base power",
 FullSimplify[(1-eps)^(1/m)==Exp[Log[1-eps]/m],m>=1&&0<eps<1]];
record["Logarithm of the base is negative",FullSimplify[Log[1-eps]<0,0<eps<1]];
record["Positive exponent times negative logarithm",
 Resolve[ForAll[{q,ell},Implies[q>0&&ell<0,q ell<0]],Reals]];
record["Exponential is at most one on negative arguments",Reduce[z<=0&&Exp[z]>1,z,Reals]===False];
am=(m-1)/(2m);
expr=1+4(u^2+am u^4)-16u^3(1+am u^2)^(3/2)+64(u^2+am u^4)^2;
s=Normal[Series[expr,{u,0,4}]];
Print["n=2 transferred polynomial: ",InputForm[Factor[s]]];
record["n=2 true third coefficient transferred",
 FullSimplify[Coefficient[s,u,4]==66-2/m,m>=1]];
record["n>=3 parameter quadratic error is within the claimed new remainder",
 FullSimplify[2>=1+3/n,n>=3]];
record["Fractional-power parameter error is absorbed in every n>=2",
 FullSimplify[2+1/n>=1+3/n,n>=2]];
record["n=2 cubic parameter error is absorbed",3>5/2];
record["Strict hierarchy of known fractional orders",
 FullSimplify[1+1/n<1+2/n<1+3/n,n>=2]];
Print["Failures: ",failures];
Exit[If[failures==0,0,1]];
