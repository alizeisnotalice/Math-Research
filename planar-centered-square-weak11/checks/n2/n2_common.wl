ClearAll[t,x,y,aa,bb,vertex,la,lb,area,rat];
rat={aa->(1-t^2)/(1+t^2),bb->2t/(1+t^2)};
vertex={0,(1-aa)/2,(1-bb)/2,(1+bb)/2,(1+aa)/2,1};
la={aa,aa,(aa+bb)/2,(aa-bb)/2,0,0};
lb={bb,bb,bb,0,0,0};
area[i_,j_]:=Together[(1+(1-vertex[[i]])(1-vertex[[j]])-la[[i]]la[[j]]-lb[[i]]lb[[j]])/.rat];
