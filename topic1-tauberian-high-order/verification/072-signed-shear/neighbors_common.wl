ClearAll[e,onegroup];
onegroup[n_,sign_,j_]:=Module[{vecs,rows,coords,inside,pred,viol,present,t0=AbsoluteTime[],cpu=TimeUsed[]},
 vecs=Select[Tuples[{-1,0,1},n],Count[#,0]==n-j&];
 rows=Table[coords=Expand[((1-sign e)IdentityMatrix[n]+sign e ConstantArray[1,{n,n}]).k];
  inside=And@@Flatten[({#>-1,#<1}& /@ coords)];
  pred=If[sign==1,Total[k]==0,Abs[Total[k]]==j&&j>=2];
  viol=Reduce[0<e<1/(4(n-1))&&If[pred,Not[inside],inside],e,Reals];
  present=Reduce[0<e<1/(4(n-1))&&inside,e,Reals];
  <|"k"->k,"predicted_inside"->pred,"actual_interior_domain"->ToString[present,InputForm],"mismatch"->!TrueQ[viol===False]|>,{k,vecs}];
 <|"n"->n,"sign"->sign,"support"->j,"kernel_id"->$KernelID,"pid"->$ProcessID,"count"->Length[rows],"predicted_count"->Count[Lookup[rows,"predicted_inside"],True],"mismatches"->Count[Lookup[rows,"mismatch"],True],"seconds"->N[AbsoluteTime[]-t0],"cpu_seconds"->N[TimeUsed[]-cpu],"max_memory_bytes"->MaxMemoryUsed[],"vectors"->rows|>];
