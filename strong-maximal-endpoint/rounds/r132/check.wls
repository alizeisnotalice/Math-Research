as=0<e<=1/4096;
tests={FullSimplify[1/(4(1-2Sqrt[2e])^2)<8/25,as],FullSimplify[(4/5)Sqrt[e](1-e)/2>e,as],FullSimplify[(2/(1-Sqrt[e])-1)^2-1==4Sqrt[e]/(1-Sqrt[e])^2,0<e<1],9^2/(10^2*2*(1-1/5))==81/160};Print[tests];If[!And@@tests,Exit[1]];
