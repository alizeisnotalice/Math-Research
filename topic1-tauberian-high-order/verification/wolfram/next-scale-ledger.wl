(* Abstract two-choice ledger only, not a realizable geometric counterexample. *)
Clear[aa,bb,v,x,xstar,value,assumptions,failures,record];
Print["VERSION: ",$Version];failures=0;
record[name_,expr_]:=Module[{r=expr},Print[name,": ",ToString[r,InputForm]];If[r=!=True,failures++]];
xstar=2 aa v^2/(1+2 bb v);value=-2 aa bb v^3/(1+2 bb v);
assumptions=aa>0&&bb>0&&0<v<1&&2 aa v^2<=1;
record["crossing_identity",FullSimplify[xstar/2-aa v^2==-bb v xstar==value,assumptions]];
record["crossing_in_domain",Resolve[ForAll[{aa,bb,v},Implies[assumptions,0<xstar<=1]],Reals]];
record["left_branch_lower_bound",Resolve[ForAll[{aa,bb,v,x},Implies[assumptions&&0<=x<=xstar,-bb v x>=value]],Reals]];
record["right_branch_lower_bound",Resolve[ForAll[{aa,bb,v,x},Implies[assumptions&&xstar<=x<=1,x/2-aa v^2>=value]],Reals]];
record["nonzero_cubic_minimax_limit",FullSimplify[Limit[value/v^3,v->0,Direction->"FromAbove"]==-2 aa bb,aa>0&&bb>0]];
Print["minimax value = ",ToString[value,InputForm]];
Print["Failures: ",failures];Exit[If[failures==0,0,1]];
