ClearAll[v,n,j];
Print["Wolfram ",$Version," PID ",$ProcessID];
rows=Flatten[Table[den=Sum[Binomial[n,k]v^k,{k,0,n-1}];
 coeffs=Normal[Series[1/den,{v,0,20}]];
 direct=Normal[Series[1/(1-(v/(1+v))^n),{v,0,n+20}]];
 Table[floor=Sum[(-1)^(j-(k-1)n)Binomial[n+j-1,j-(k-1)n],{k,1,1+Floor[j/n]}];
 <|"n"->n,"j"->j,"coefficient"->Coefficient[coeffs,v,j],"pass"->TrueQ[Coefficient[coeffs,v,j]==floor==Coefficient[direct,v,n+j]]|>,{j,0,20}],{n,2,6}],1];
Print["cases: ",Length[rows]," all pass: ",And@@Lookup[rows,"pass"]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"ell_wl.json"}],rows,"RawJSON"];
