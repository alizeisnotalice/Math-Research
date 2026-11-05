ClearAll[coefc,ap,aq,onecase,w,z];
coefc[q_,i_]:=If[i==0,1,Product[t-q,{t,1,i}]/(i+1)!];
ap[p_,kk_,q_]:=Module[{arr=ConstantArray[0,kk+1]},arr[[1]]=1;Do[arr[[k+1]]=FullSimplify[Sum[((p+1)i-k)coefc[q,i]arr[[k-i+1]],{i,1,k}]/k],{k,1,kk}];arr[[kk+1]]];
aq[q_,z_]:=Which[q==0,-Log[1-z]/z,q==1,1,True,(1-(1-z)^q)/(q z)];
onecase[n_,r_,q_]:=Module[{A,mat,direct,expected,mis={},j,J,k,t0=AbsoluteTime[]},
 A=Normal[Series[aq[q,z],{z,0,Floor[r/n]}]];
 Do[direct=Normal[Series[w^(n+j)(A/.z->w^n)^((n+j)/n),{w,0,n+r}]];
  Do[k=(J-j)/n;expected=If[J>=j&&IntegerQ[k],ap[(n+j)/n,k,q],0];
   If[!TrueQ[FullSimplify[Coefficient[direct,w,n+J]-expected]==0],AppendTo[mis,{j,J}]],{J,0,r}],{j,0,r}];
 <|"n"->n,"r"->r,"q"->ToString[q,InputForm],"kernel_id"->$KernelID,"pid"->$ProcessID,"seconds"->N[AbsoluteTime[]-t0],"basis_assertions"->(r+1)^2,"mismatches"->Length[mis],"pass"->(mis=={})|>];
