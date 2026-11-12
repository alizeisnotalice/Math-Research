Clear[n,w,l,q];
checks={Resolve[ForAll[{n,w,l,q},Implies[n>=1&&w>0&&0<q<1&&l>10n w/q,2n w/(l-w)<q&&l>w+2w/q]],Reals]};
Do[With[{dirs=Subsets[Range[n],{d}],mapped=Sort[Take[#,d]]& /@ Permutations[Range[n]]},AppendTo[checks,Sort[Keys[Counts[mapped]]]===Sort[dirs]&&Union[Values[Counts[mapped]]]==={Factorial[d]Factorial[n-d]}]],{n,2,6},{d,1,n-1}];
Print[{Length[checks],And@@checks}];If[!And@@checks,Exit[1]];
