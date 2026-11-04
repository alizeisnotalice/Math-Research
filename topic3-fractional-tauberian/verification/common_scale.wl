ClearAll[a,t,l];hh=(1+2a^2)/Sqrt[1+4a^2]-a;
inv=InverseSeries[Series[1-hh,{a,0,10}]];
Print[<|"inverse_series"->ToString[inv/.a->t,InputForm],"H_series"->ToString[Series[hh,{a,0,10}],InputForm],"derivative"->FullSimplify[D[hh,a]],"layer_identity"->FullSimplify[(3+(Sqrt[1+4a^2]-2a)^2)/(4Sqrt[1+4a^2])-hh],"version"->$Version|>];
