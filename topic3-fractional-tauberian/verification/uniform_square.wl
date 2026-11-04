(* Exact identities for the all-order envelope and the side-window coefficients. *)
ClearAll[b,t,u,j,l,q];
r={ (1-(1-t)^(1/b))/2,
 1-t u-(1-t)(1-2t u)^(1-b),
 1+t u-(1-t)(1+2t u)^(2-b),
 ((1-t)^(-1/(2-b))-1)/2};
a[j_]:={(-1)^(j+1) Binomial[1/b,j]/2,
 -Binomial[1-b,j](-2u)^j+Binomial[1-b,j-1](-2u)^(j-1)-u KroneckerDelta[j,1],
 -Binomial[2-b,j](2u)^j+Binomial[2-b,j-1](2u)^(j-1)+u KroneckerDelta[j,1],
 (-1)^j Binomial[-1/(2-b),j]/2};
res=Table[FullSimplify[FunctionExpand[Coefficient[Normal[Series[r,{t,0,6}]],t,j]-a[j]],
 Assumptions->0<b<2],{j,1,6}];
low=Normal[Series[r[[{1,2}]]/.b->1/2+l t,{t,0,2}]];
high=Normal[Series[r[[{4,3}]]/.b->3/2+l t,{t,0,2}]];
resLow=Simplify[low-{t+(-1/2-2l)t^2,t+(-(1+2l)u+u^2/2)t^2}];
resHigh=Simplify[high-{t+(3/2+2l)t^2,t+((1+2l)u+u^2/2)t^2}];
resG={Simplify[4 Integrate[(1-u)(q-u),{u,0,q}]-(2q^2-2q^3/3)],
 Simplify[4 Integrate[(1-u)(q-u),{u,0,1}]-(2q-2/3)]};
(* q=4 gives an independent exact omitted-coefficient check for N=1..6. *)
rd=((1-t)^(-4)-1)/2;
sharp=Table[Simplify[Coefficient[Normal[Series[4(rd^2-(Normal[Series[rd,{t,0,k}]])^2),{t,0,k+2}]],t,k+2]
 -8 Pochhammer[4,k+1]/Factorial[k+1]],{k,1,6}];
If[sharp=!=ConstantArray[0,6],Exit[1]];
result=<|"kind"->"exact-symbolic-identities","version"->$Version,
 "coefficient_residuals_N6"->res,"low_window_residuals"->resLow,
 "high_window_residuals"->resHigh,"G_integral_residuals"->resG,"sharpness_omitted_coefficient_residuals"->sharp|>;
Print[result];
If[Flatten[res]=!=ConstantArray[0,24]||resLow=!={0,0}||resHigh=!={0,0}||resG=!={0,0},Exit[1]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"..","results","uniform_square_wolfram.json"}],result,"RawJSON"];
