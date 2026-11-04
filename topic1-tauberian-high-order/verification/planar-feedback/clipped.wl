ClearAll[eta,t,x,y,h];
Print["Wolfram version: ", $Version];
Print["Domain: 0<eta<1/2; target arms x<=eta(1-y) or y<=eta(1-x), unit square."];
checks={
 Reduce[0<eta<1/2 && 0<=t<=1-eta && eta (1-t)<eta^2,{eta,t},Reals],
 Reduce[0<eta<1/2 && 1-eta<t<=1 && t<=eta,{eta,t},Reals],
 Reduce[0<eta<1/2 && 0<=t<=eta^2 && 1-t/eta<eta(1-t),{eta,t},Reals]};
Print["branch violations (all False): ",checks];
first=FullSimplify[(Integrate[(1-t)eta^2,{t,0,1-eta}]+Integrate[(1-t)eta(1-t),{t,1-eta,1}])/eta,0<eta<1/2];
second=FullSimplify[Integrate[(eta^2-t)(1-t/eta),{t,0,eta^2}]/eta,0<eta<1/2];
formal=1/3+3 eta^2/2-eta^6/6;
Print["first true/eta: ",first];Print["second true/eta: ",second];
Print["true total/eta: ",Expand[first+second]];Print["formal total: ",formal];
Print["gap nonpositive counterexample: ",Reduce[0<eta<1/2 && formal-first-second<=0,eta,Reals]];
Print["limits {true,formal}: ",{Limit[first+second,eta->0,Direction->"FromAbove"],Limit[formal,eta->0,Direction->"FromAbove"]}];
(* Unit symmetric weighted union: 2*one arm minus exact overlap. *)
arm=Integrate[Integrate[2-x-y,{x,0,h(1-y)}],{y,0,1},Assumptions->0<h<1/4];
overlap=Integrate[Integrate[2-x-y,{y,0,h(1-x)}],{x,0,h/(1+h)},Assumptions->0<h<1/4]+Integrate[Integrate[2-x-y,{y,0,1-x/h}],{x,h/(1+h),h},Assumptions->0<h<1/4];
weighted=FullSimplify[2 arm-overlap,0<h<1/4];
Print["symmetric weighted union: ",weighted];
Print["difference from 5h/3: ",Factor[weighted-5h/3]];
Print["positive difference counterexample: ",Reduce[0<h<=1/4 && weighted>5h/3,h,Reals]];
Print["quadratic correction limit: ",Limit[(weighted-5h/3)/h^2,h->0,Direction->"FromAbove"]];
Print["checks complete"];
