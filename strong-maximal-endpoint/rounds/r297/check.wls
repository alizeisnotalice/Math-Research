Clear[price];price[d_,t_,u_]:=Sum[(t+u)^j/j!,{j,0,d-1}]/Sum[t^j/j!,{j,0,d-1}];
checks=Flatten[Table[With[{v=price[d,t,u],upper=(1+u/t)^(d-1)},upper/2<=v<=upper],{d,{2,4,8,16,32}},{t,{2d,d^3}},{u,{0,1,d,d^2}}]];
Print[{Length[checks],And@@checks}];If[!And@@checks,Exit[1]];
