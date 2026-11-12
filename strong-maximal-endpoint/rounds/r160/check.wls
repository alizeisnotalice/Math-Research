a={{1,1},{0,1}};
k=Sum[a[[i,j]]a[[l,h]](1-a[[i,h]])(1-a[[l,j]]),{i,2},{j,2},{l,2},{h,2}]/16;
tests={k==0,Total[Flatten[a]]/4==3/4,2(1/4)^2==1/8};Print[tests];If[!TrueQ[And@@tests],Exit[1]];
