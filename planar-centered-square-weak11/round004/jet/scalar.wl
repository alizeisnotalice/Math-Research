Print["PID ",$ProcessID," version ",$Version];
base=If[Length[$ScriptCommandLine]>1,Last[$ScriptCommandLine],FileNameJoin[{DirectoryName[$InputFileName],"generated"}]];
z=Import[FileNameJoin[{base,"jet.json"}],"RawJSON"];
d=ToExpression[z["score_dA_dt"]];gap=ToExpression[z["minimum_strict_mass_gap"]];until=ToExpression[z["positive_radicand_until"]];
Print["score positive ",d>0];Print["uniform strict gap positive ",gap>0];Print["t=1e-9 radicands positive ",1/10^9<until];
Print["calibration preserves strict gap: ",FullSimplify[(M+t w)-4(r^2+t w/4)==M-4r^2]];
Print["phys derivative scaling: ",FullSimplify[(D/(8L))==((L D/8)/L^2),Assumptions->L>0]];
