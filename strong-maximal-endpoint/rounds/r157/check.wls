cross[a_]:=With[{n=Length[a]},Sum[a[[i,j]]a[[k,l]](1-a[[k,j]])(1-a[[i,l]]),{i,n},{k,n},{j,n},{l,n}]/n^4];
tests={cross[IdentityMatrix[2]]==1/8,cross[{{1,1},{0,1}}]==0};Print[tests];If[!TrueQ[And@@tests],Exit[1]];
