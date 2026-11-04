(* Multi-depth scalar ledger. No geometric realization is asserted here. *)
Clear[k,c,b,v,A,f,weights,zstar,fstar,failures,record];
Print["VERSION: ",$Version];failures=0;
record[name_,expr_]:=Module[{r=expr},Print[name,": ",ToString[r,InputForm]];If[r=!=True,failures++]];
Do[
 weights=Table[(b v/c)^(k-j),{j,1,k}];
 record["telescoping_k_"<>ToString[k],
  Simplify[Sum[weights[[j]](c f[j]-b v f[j-1]),{j,1,k}]==c f[k]-b v (b v/c)^(k-1) f[0],c>0&&b>0&&v>0]];
 zstar=-b v (b v/c)^(k-1) A v/Total[weights];
 fstar[0]=A v;
 Do[fstar[j]=Together[(zstar+b v fstar[j-1])/c],{j,1,k}];
 record["attainment_terminal_k_"<>ToString[k],Simplify[fstar[k]==0,c>0&&b>0&&v>0&&A>0]];
 record["nonzero_order_k_"<>ToString[k],
  Simplify[Limit[zstar/v^(k+1),v->0,Direction->"FromAbove"]==-A b^k/c^(k-1),c>0&&b>0&&A>0]],
 {k,1,6}];
Print["Failures: ",failures];Exit[If[failures==0,0,1]];
