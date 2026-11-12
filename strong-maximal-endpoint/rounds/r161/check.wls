a={{1,1},{1,0},{0,1}};k=Sum[a[[i,j]]a[[l,h]](1-a[[i,h]])(1-a[[l,j]]),{i,3},{j,2},{l,3},{h,2}]/36;
tests={k==1/18,Total[Flatten[a]]/6-1/2==1/6};Print[tests];If[!TrueQ[And@@tests],Exit[1]];
