p=1/2;
vals=Table[Module[{d,m,balance},d[t_]:=Boole[t<1]/2+Boole[t<h]/4;
m=Integrate[d[t]-p,{t,0,1}];balance=Integrate[p-d[t],{t,1,3}];
{m,balance,2m==balance}],{h,{3,4}}];
Print[vals];If[!And@@vals[[All,3]],Exit[1]];
