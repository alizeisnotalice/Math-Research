ClearAll[a,b,c,u,v,w,h];
sol=Maximize[{h,0<=a<=3/10&&0<=b<=3/10&&0<=c<=3/10&&0<=u<=1/10&&0<=v<=1/10&&0<=w<=1/10&&a+b+c+u+v+w<=1&&a+b+u>=16/25+h&&a+c+v>=16/25+h&&b+c+w>=16/25+h},{a,b,c,u,v,w,h},Reals];
Print[<|"exact_primal_optimum"->First[sol],"dual_value"->(1/3+9/10*(2/3-1/3)-16/25),"individual_margin"->(7/10-16/25),"pair_margin"->(1+3/10-2*16/25),"version"->$Version|>];
