over[g_,l_,r_]:=Max[0,Min[g,l]+Min[g,r]-g];
g={1/4,3/2};
actual=over[g[[1]],1,1]+over[g[[2]],7/4,1];
local=over[g[[1]],1,1]+over[g[[2]],1,1];
margin=(15/4)(3)/(33/8)^2-1/2;
Print[{actual,local,actual-local,margin}];
If[actual!=5/4||local!=3/4||margin<=0,Exit[1]];
