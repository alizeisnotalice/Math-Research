(* Exact coefficient identities; bounded-parameter remainder is proved in the paper. *)
ClearAll[t,u,l,m,a];
outerB=Normal[Series[1-2t u-(1-t)(1-2t u)^(-l t),{t,0,3}]];
outerC=Normal[Series[1-(1-t)(1+2t u)^(1-l t),{t,0,3}]];
expectedB=t(1-2u)-2l u t^2+2l(u-u^2)t^3;
expectedC=t(1-2u)+2(1+l)u t^2-2l(u-u^2)t^3;
resOuter=Simplify[{outerB-expectedB,outerC-expectedC}];
innerB=Normal[Series[1-2t u-(1-t)(1-2t u)^(t/2-m t^2),{t,0,3}]];
innerC=Normal[Series[1-(1-t)(1+2t u)^(1+t/2-m t^2),{t,0,3}]];
resInner=Simplify[{innerB-(t(1-2u)+u t^2+u(u-1-2m)t^3),
 innerC-(t(1-2u)+u t^2-u(u-1-2m)t^3)}];
(* For a in [0,1/2], split the absolute-value integral at u=a. *)
q=1/2+8(Integrate[u(a-u),{u,0,a}]+Integrate[u(u-a),{u,a,1/2}]);
resIntegral=Simplify[q-(5/6-a+8a^3/3)];
result=<|"kind"->"exact-symbolic-identities","version"->$Version,
 "outer_residuals"->resOuter,"inner_residuals"->resInner,
 "split_integral_residual"->resIntegral,"middle_Q"->ToString[Expand[q],InputForm],
 "assumptions"->"real bounded lambda/mu, t tends to zero from above; middle integral 0<=a<=1/2"|>;
Print[result];
If[resOuter=!={0,0}||resInner=!={0,0}||resIntegral=!=0,Exit[1]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"..","results","critical_square_wolfram.json"}],result,"RawJSON"];
