ClearAll["Global`*"];
Print["version=",$Version]; Print["pid=",$ProcessID];
Print["measurePrimal=",Minimize[{m1+m2+m12,m1>=0&&m2>=0&&m12>=0&&m1+m12>=1&&m2+m12>=1},{m1,m2,m12}]];
Print["measureDual=",Maximize[{z1+z2,z1>=0&&z2>=0&&z1+z2<=1},{z1,z2}]];
Print["densityPrimal=",Minimize[{m1+m2,m1>=1&&m2>=1},{m1,m2}]];
Print["densityDual=",Maximize[{z1+z2,0<=z1<=1&&0<=z2<=1},{z1,z2}]];
eta=q/(p+q);
Print["sqrtContinuityOptimizerResidual=",FullSimplify[p^2/(1-eta)+q^2/eta-(p+q)^2,p>0&&q>0]];
