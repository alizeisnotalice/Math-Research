(* For equal balanced blocks at the median f(v)=1-v. *)
vals=Table[Total[Table[Integrate[Abs[v-(2i-1)/(2k)],{v,(i-1)/k,i/k}],{i,k}]],{k,{2,4,8}}];
Print[vals];If[!TrueQ[vals=={1/8,1/16,1/32}],Exit[1]];
