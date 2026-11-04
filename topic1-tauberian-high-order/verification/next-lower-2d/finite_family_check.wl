SetDirectory[DirectoryName[$InputFileName]];Print["VERSION ",$Version," PID ",$ProcessID];
Clear[decExact];decExact[s_String]:=Module[{p=StringSplit[s,RegularExpression["[eE]"]],b,x=0,sgn=1,dp},b=First[p];If[Length[p]>1,x=ToExpression[p[[2]]]];If[StringStartsQ[b,"-"],sgn=-1;b=StringDrop[b,1]];dp=If[StringContainsQ[b,"."],StringLength[Last[StringSplit[b,"."]]],0];sgn FromDigits[StringReplace[b,"."->""]]10^(x-dp)];
r=Import["finite_family.json","RawJSON"];u=ToExpression[r["u"]];kap=ToExpression[r["kappa"]];lam=ToExpression[r["lambda_rational"]];w=u/(1+2u)+lam u^6+u^8;rho=1-2w;e=kap u^3;a=4u^2;det=1-e^2;m0=(2u/(1+2u))^2;mass=a rho^2;dv=4w^2-mass;
h=hs/.FindRoot[a(2rho hs+hs^2)==e^2,{hs,N[e^2/(2a rho),210]},WorkingPrecision->210,AccuracyGoal->170,PrecisionGoal->170];
area=2h NIntegrate[e-a(2rho(h z)+(h z)^2)/e,{z,0,1},WorkingPrecision->170,AccuracyGoal->130,PrecisionGoal->130];
budget=dv-2area;G=(mass-det m0)/((det-mass)(1-m0)a(2u)^6);Nn=2^90;fd=((2Nn-5)/(2Nn))^2 mass/(det-mass)-m0/(1-m0);radius=rho+w;margin=u^2-(mass+2e^2)/(2radius)^2;
idx=FromDigits[StringDrop[r["hstar"]["lower_hex"],2],16];inside[root_,bounds_]:=TrueQ[decExact[bounds[[1]]]<=root<=decExact[bounds[[2]]]];
checks={"hBracket"->TrueQ[(w-e)idx/2^400<=h<=(w-e)(idx+1)/2^400],"area"->inside[area,r["A_bounds"]],"budget"->inside[budget,r["D_minus_2A"]],"G"->inside[G,r["G"]],"strictBudget"->TrueQ[budget>0],"centerMargin"->TrueQ[margin==ToExpression[r["center"]["strict_margin"]]&&margin>0],"finiteGap"->TrueQ[fd==ToExpression[r["finite_halo_minus_Ln"]["exact_fraction"]]&&fd>0]};
Print["CHECKS ",InputForm[checks]];Print["NUMERIC ",InputForm[{"G"->N[G,30],"budget"->N[budget,30],"h"->N[h,30]}]];
