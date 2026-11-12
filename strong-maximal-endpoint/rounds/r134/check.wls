tests={FullSimplify[1/(1-4h)+2k h<=1+1/(10k),k>=2&&0<h<=1/(100k^2)],E^(1/10)<5/4,15/8/2<1};
vals=Table[N[(k!)^(1/k)/k ((3/2)^(1/k)-(5/4)^(1/k)) k,16],{k,{2,4,8,16,32,64,128,256}}];
Print[{tests,vals}];If[!And@@tests,Exit[1]];
