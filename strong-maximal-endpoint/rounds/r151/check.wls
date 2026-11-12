e=1/32;t=1/2;
row[u_]:={{0,u-e},{u,u+e}};
len[r_]:=Total[(Last[#]-First[#])& /@ r];
over[i_,j_]:=Max[0,Min[i[[2]],j[[2]]]-Max[i[[1]],j[[1]]]];
c=Total[Flatten[Table[over[i,j],{i,row[t]},{j,row[t+e]}]]];
vals={len[row[t]]-c,len[row[t+e]]-c,e(1/2-e)/2,1/2+e^2};
Print[vals];If[!TrueQ[Take[vals,2]=={e,2e}],Exit[1]];
