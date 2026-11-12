tests={FullSimplify[2/(2-x)>8/15&&1/(1-x)<=2/5&&1/(2-x)<=2/7,-7/4<x< -3/2],(7/4-3/2)(1-1/2)/2==1/16};Print[tests];If[!TrueQ[And@@tests],Exit[1]];
