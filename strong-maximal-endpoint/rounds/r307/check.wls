Clear[x];checks={Resolve[ForAll[x,Implies[x>=3,Equivalent[Max[1/(x-2),2/x]>1/2,x<4]]],Reals]};
Do[With[{p=Permutations[Range[n]]},AppendTo[checks,Union[Values[Counts[First/@p]]]==={Factorial[n-1]}]],{n,2,7}];
Print[{Length[checks],And@@checks}];If[!And@@checks,Exit[1]];
