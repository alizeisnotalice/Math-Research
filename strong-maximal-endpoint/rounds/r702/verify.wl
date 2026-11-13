poly=Exp[-x]Sum[x^j/j!,{j,0,63}];
out=<|"gamma_tail_derivative"->ToString[FullSimplify[D[poly,x]+Exp[-x]x^63/63!],InputForm],
"gamma_tail_at_zero"->ToString[poly/.x->0,InputForm],
"gamma_tail_at_infinity"->ToString[Limit[poly,x->Infinity],InputForm],
"modular_derivative"->ToString[Factor[D[h(1+Log[h])^63,h]],InputForm]|>;
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
