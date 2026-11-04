ClearAll[u,a,b,c,d,e,f,g,j];
Print["Mathematica: ",$Version];
mat={{1,u^3 a,u^3 b},{u^3 c,1,u^3 d},{u e+u^3 f,u j+u^3 g,1}};
(* General order-three last-row perturbations use independent variable g. *)

pair=mat[[1,2]]mat[[2,1]]+mat[[1,3]]mat[[3,1]]+mat[[2,3]]mat[[3,2]];
cycles=mat[[1,2]]mat[[2,3]]mat[[3,1]]+mat[[1,3]]mat[[2,1]]mat[[3,2]];
Print["exact det identity residual: ",Expand[Det[mat]-(1-pair+cycles)]];
Print["cycle polynomial: ",Factor[cycles]];
Print["coefficients powers 0..6: ",Table[Coefficient[cycles,u,j],{j,0,6}]];


Print["explicit three-axis box intersection: ",Expand[(2 b u^3)^2 (2 s u)]];
ClearAll[s];
tmat={{1,-u^3,-u^3},{-u^3,1,-u^3},{-u,-u,1}};
shift=Expand[tmat.{1,1,1}];
Print["explicit affine shift b=s=1: ",shift];
Print["determinant: ",Factor[Det[tmat]]];
Print["guard violation 0<u<1/8: ",Reduce[0<u<1/8 && (Det[tmat]<=0 || 2u^3>=1 || 2u>=1),u,Reals]];
Print["axis overlap product: ",Expand[Times@@(1-shift)]];

Print["complete"];
