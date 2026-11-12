(* Equal blocks and median: f(u)=1-u gives exact loss from the triangles. *)
checks=Table[With[{v=Range[0,m+1]/(m+1)},1/2-Total[Differences[v] Most[v]]==1/(2(m+1))],{m,{1,2,4,8,16,32,64,128}}];Print[checks];If[!TrueQ[And@@checks],Exit[1]];
