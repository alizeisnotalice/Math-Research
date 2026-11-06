(* r044 exact k=2 supplement for the fixed r043 tensor dictionary. *)
root442=DirectoryName[$InputFileName];out442=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{$TemporaryDirectory,"tensor_k2_exact44.json"}]];
c442=69/100;
(* m_A is the best V/(|E|A) among r043 rows of amplitude A. *)
types442={{2,11/8,c442},{3,7/3,3c442/2},{4,3,2c442},{8,4,3c442},{16,15/4,4c442}};
pairs442=Flatten[Table[Module[{x=types442[[i]],y=types442[[j]],ub},ub=x[[2]]y[[2]]/(1+x[[3]]+y[[3]])^3;
 <|"amplitudes"->{x[[1]],y[[1]]},"m_product"->ToString[x[[2]]y[[2]],InputForm],"rational_log_lower_sum"->ToString[x[[3]]+y[[3]],InputForm],"rational_upper_bound"->ToString[ub,InputForm],"upper_decimal"->N[ub,30],"strictly_less_than_1"->(ub<1)|>],{i,Length[types442]},{j,i,Length[types442]}],1];
best442=First@SortBy[pairs442,-ToExpression[#"rational_upper_bound"]&];
If[!And@@pairs442[[All,"strictly_less_than_1"]],Print["k2 certificate failed"];Exit[1]];
res442=<|"status"->"r044 exact rational k=2 tensor-pair certificate; fixed r043 dictionary only",
 "m_by_amplitude"-><|"2"->"11/8","3"->"7/3","4"->"3","8"->"4","16"->"15/4"|>,
 "log_lower_bounds"-><|"Log[2]"->"69/100","Log[3]"->"207/200","Log[4]"->"69/50","Log[8]"->"207/100","Log[16]"->"69/25"|>,
 "justification"->"Log[2]>69/100 follows from the r044 rational Taylor certificate. Log[3]>(3/2)Log[2] because 9>8; powers of 2 give the remaining bounds.",
 "k2_formula"->"L_4=m_i m_j/(1+Log[A_i]+Log[A_j])^3 <= listed rational_upper_bound",
 "all_15_unordered_pairs"->pairs442,"largest_rational_upper_bound"->best442,
 "consequence"->"Every two-factor heterogeneous product from the fixed 20 witnesses is strictly below 1. With r044's B_k proof for k>=3, every fixed-dictionary product witness with k>=2 is strictly below 1."|>;
Export[out442,res442,"RawJSON"];Print[<|"largest"->best442,"all_pass"->And@@pairs442[[All,"strictly_less_than_1"]]|>];
