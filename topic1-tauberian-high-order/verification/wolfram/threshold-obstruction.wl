(* Exact scalar checks for the finite-tail obstruction, with no run-time limits. *)
Clear[J,k,m,v,dd,n];
Print["Wolfram version: ", $Version];
ass=Element[J,Integers] && J>=2;
wsum=FullSimplify[Sum[1/(J 2^k),{k,1,J}],Assumptions->ass];
Print["wsum = ",InputForm[wsum]];
Print["mean = ",InputForm[FullSimplify[1-wsum+Sum[2^k/(J 2^k),{k,1,J}],Assumptions->ass]]];
tail=FullSimplify[Sum[1/(J 2^k),{k,m+1,J}],Assumptions->ass && Element[m,Integers] && 0<=m<=J-1];
Print["tail = ",InputForm[tail]];
Print["LF lower = ",InputForm[FullSimplify[2^m tail,Assumptions->ass && Element[m,Integers] && 0<=m<=J-1]]];
Print["budgetMargin = ",InputForm[FullSimplify[(2^(n+2)dd)(1-dd)/2^n-dd,Assumptions->n>=2 && dd>0]]];
Print["costRatioAtL1 = ",InputForm[FullSimplify[(v^2+2v(1-v)/J)/(v/J)]]];
