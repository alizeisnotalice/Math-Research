(* Exact local algebra and simplex moments; the uniform theorem is proved analytically. *)
ClearAll[r,t,b,s,m,k,n,u];
res={};
Do[
 z=Prepend[Array[u,dim-1],1];
 Do[
  fm=(1-2r z[[i]])^(b-i+1) Times@@Table[1-r(z[[i]]+z[[j]]),{j,1,i-1}];
  fp=(1+2r z[[i]])^(b-dim) Times@@Table[1-r(z[[j]]-z[[i]]),{j,1,i-1}];
  amin=Sum[z[[j]],{j,1,i-1}]+(2b-i+1)z[[i]];
  aplus=Sum[z[[j]],{j,1,i-1}]+(2(dim-b)-i+1)z[[i]];
  ys=Table[z[[j]]-z[[i]],{j,1,i-1}];
  bm=-Total[ys^2]/2-2z[[i]] Total[ys]-2b z[[i]]^2;
  bp=-Total[ys^2]/2+2(dim-b)z[[i]]^2;
  AppendTo[res,Simplify[Coefficient[Normal[Series[Log[fm],{r,0,2}]],r,1]+amin]];
  AppendTo[res,Simplify[Coefficient[Normal[Series[Log[fm],{r,0,2}]],r,2]-bm]];
  AppendTo[res,Simplify[Coefficient[Normal[Series[Log[fp],{r,0,2}]],r,1]+aplus]];
  AppendTo[res,Simplify[Coefficient[Normal[Series[Log[fp],{r,0,2}]],r,2]-bp]];
 ,{i,1,dim}];
,{dim,1,5}];
(* Beta identities with m=n-k, including k=0 by continuity. *)
moments={Beta[m+2,k+1]/(Beta[m,k+1]s^2)-m(m+1)/(s^2(k+m+1)(k+m+2)),
 k Beta[m+1,k+2]/((k+1)s Beta[m,k+1])-k m/(s(k+m+1)(k+m+2)),
 2k Beta[m,k+3]/((k+1)(k+2)Beta[m,k+1])-2k/((k+m+1)(k+m+2))};
momres=FullSimplify[FunctionExpand[moments],Assumptions->m>0&&k>=0&&s>0];
minus[n_,k_,s_]:=n/2-(k+(n-k)(n+k+1)/s)/(n+1);
plus[n_,k_,s_]:=n/2+((n-k)(n-k+1)/s-k)/(n+1);
jump=Simplify[minus[n,k,k]-minus[n,k-1,k]-1/(n+1)];
half=Simplify[plus[n,n-1,n]-minus[n,n-1,n]-2/n];
coeff[n_,b_]:=Module[{c=Min[b,n-b],kk,ll,tt},
 kk=If[b<n/2,Floor[2b],Ceiling[2(n-b)]-1];
 ll=2^n/(kk!(2c)^(n-kk));
 tt=ll If[b<n/2,minus[n,kk,2c],plus[n,kk,2c]];
 <|"n"->n,"beta"->ToString[b,InputForm],"k"->kk,"L"->ToString[ll,InputForm],"T"->ToString[tt,InputForm],"L_numeric"->N[ll,30],"T_numeric"->N[tt,30]|>];
rows=Join[Table[coeff[1,j/4],{j,1,3}],Table[coeff[2,j/4],{j,1,7}],
 Table[coeff[3,j/4],{j,1,11}],Table[coeff[4,j/2],{j,1,7}],Table[coeff[5,j/2],{j,1,9}]];
(* Finite exact Lagrange inversions through order six at rational directions. *)
inverseRes={};
Do[
 zz=Table[(dim-j+1)/dim,{j,1,dim}];
 Do[
  Do[
   ff=If[sgn<0,(1-2r zz[[i]])^(bb-i+1) Times@@Table[1-r(zz[[i]]+zz[[j]]),{j,1,i-1}],
     (1+2r zz[[i]])^(bb-dim) Times@@Table[1-r(zz[[j]]-zz[[i]]),{j,1,i-1}]];
   cc=Table[Coefficient[Normal[Series[(r/(1-ff))^jj,{r,0,jj-1}]],r,jj-1]/jj,{jj,1,6}];
   inv=Sum[cc[[jj]]t^jj,{jj,1,6}];
   AppendTo[inverseRes,Normal[Series[(ff/.r->inv)-(1-t),{t,0,6}]]//Simplify];
  ,{i,1,dim},{sgn,{-1,1}}];
 ,{bb,{(dim-1)/2,dim/2}}];
,{dim,{3,4,5}}];
sharpRes=Flatten[Table[With[{rr=((1-t)^(-4)-1)/2},
 Simplify[Coefficient[Normal[Series[(2rr)^dim-(2Normal[Series[rr,{t,0,nn}]])^dim,{t,0,dim+nn}]],t,dim+nn]
 -dim 4^(dim-1)Pochhammer[4,nn+1]/Factorial[nn+1]]],{dim,1,5},{nn,1,6}]];
If[!And@@Thread[inverseRes==0]||!And@@Thread[sharpRes==0],Exit[1]];
result=<|"kind"->"exact-local-identities-and-coefficient-table","version"->$Version,
 "contact_expansion_residuals_count"->Length[res],"contact_expansion_residuals_all_zero"->And@@Thread[res==0],
 "beta_moment_residuals"->momres,"inner_jump_residual"->jump,"half_jump_residual"->half,"Lagrange_order6_cases"->Length[inverseRes],"Lagrange_residuals_all_zero"->And@@Thread[inverseRes==0],"sharp_truncation_coefficient_checks"->Length[sharpRes],"sharp_truncation_residuals_all_zero"->And@@Thread[sharpRes==0],"rows"->rows|>;
Print[KeyDrop[result,{"rows"}]];
If[!And@@Thread[res==0]||momres=!={0,0,0}||jump=!=0||half=!=0,Exit[1]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"..","results","higher_dimensional_wolfram.json"}],result,"RawJSON"];
