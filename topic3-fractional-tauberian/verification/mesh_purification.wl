ClearAll[h,u,i,j,nn];
coeffChecks=Table[And@@Map[#>=0&,Flatten[CoefficientList[Expand[4n h(u+4h)^(n-1)-((u+4h)^n-u^n)],{h,u}]]],{n,1,8}];
mass=FullSimplify[Sum[(2i+1)(2j+1)/(4nn^4),{i,0,nn-1},{j,0,nn-1}],Assumptions->Element[nn,Integers]&&nn>=1];
Print[<|"boundary_polynomial_nonnegative_n1_to8"->coeffChecks,"density_mesh_mass"->mass,"finite_margin"->FullSimplify[(1-(1/2-1/nn)-1/2)/2],"version"->$Version|>];
