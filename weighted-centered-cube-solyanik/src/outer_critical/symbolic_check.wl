Print["Wolfram version: ",$Version];Clear[tau,u,R];
pair[tau_] := 3/(4 tau)*(5/3)^tau;
checks=FullSimplify[{pair[1/2]==Sqrt[15]/2,
 pair[1/4]^4-(17/5)^4==854/625,
 pair[1/2]^2-(242/125)^2==119/62500}];
Print["Pair certificates: ",checks];If[checks=!={True,True,True},Exit[1]];
ustar=(Sqrt[2]-1)^2;
half=FullSimplify[{Sqrt[ustar]+ustar/2==1/2,
 (1+Sqrt[ustar])/((1+ustar)/2)==1+Sqrt[2]}];
Print["Exact A1 comparison tau=1/2: ",half];
If[half=!={True,True},Exit[1]];
Do[uc=u/.FindRoot[u^(1-ex)+(1-ex)u==ex,{u,1/10,0,1},
 WorkingPrecision->60];
 Print["tau=",ex,", pair lower=",N[pair[ex],30],
 ", comparison upper=",N[2 uc^(ex-1),30]],{ex,{1/4,1/2,3/4,9/10}}];
Print["m=2 outer pole lower tau=1/4: ",N[pair[1/4]/2,30]];
