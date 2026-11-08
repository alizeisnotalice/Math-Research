(* Compare independently evaluated conditional Fourier and finite-mixture formulas. *)
args=Rest[$ScriptCommandLine];If[Length[args]<3,Print["Need conditional.json mixture.json output.json"];Exit[2]];
f=Import[args[[1]],"RawJSON"];g=Import[args[[2]],"RawJSON"];
key66[r_]:=ToString[r["b"]]<>"|"<>ToString[r["kind"]]<>"|"<>r["kappa"];
fg=GroupBy[f["rows"],key66];gg=GroupBy[g["rows"],key66];keys=Intersection[Keys[fg],Keys[gg]];
rows=Table[fr=Last[SortBy[fg[k],#["refined"]&]];gr=Last[SortBy[gg[k],#["precision"]&]];
<|"key"->k,"CF_absolute_difference"->ToString[N[Abs[ToExpression[fr["CF"]]-ToExpression[gr["CF"]]],30],InputForm],"derivative_absolute_difference"->ToString[N[Abs[ToExpression[fr["derivative"]]-ToExpression[gr["derivative"]]],30],InputForm]|>,{k,keys}];
cfmax=Max[ToExpression[#["CF_absolute_difference"]]&/@rows];dmax=Max[ToExpression[#["derivative_absolute_difference"]]&/@rows];
Export[args[[3]],<|"case_count"->Length[rows],"max_CF_absolute_difference"->ToString[cfmax,InputForm],"max_derivative_absolute_difference"->ToString[dmax,InputForm],"rows"->rows,"scope"->"Independent formulas agree numerically; no directed total-error certificate"|>,"RawJSON"];Print["Compared ",Length[rows]," independent conditional cases"];
If[Length[rows]<1||cfmax>10^-25||dmax>10^-25,Exit[3]];
