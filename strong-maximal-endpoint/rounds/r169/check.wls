(* At the median the exact covariance is ArcSin[z]/(2 Pi). *)
he[j_]:=2^(-j/2) HermiteH[j,0];
checks=Table[FullSimplify[SeriesCoefficient[ArcSin[z]/(2Pi),{z,0,j}]==he[j-1]^2/(2Pi j!)],{j,1,12}];
Print[{Length[checks],And@@checks}];If[!TrueQ[And@@checks],Exit[1]];
