h=1/4;d=1/32;p=1/2+d-d^2/2;
shared=Integrate[Min[1,1+d-v],{v,0,h}];
tests={shared==h-(h-d)^2/2,p>1/2,(3/8-1/8)^2==1/16,2-h==7/4};Print[{tests,2p-shared,p-1/2}];If[!TrueQ[And@@tests],Exit[1]];
