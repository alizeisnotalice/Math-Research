s[h_]:=Sum[(4/5)^j,{j,1,h}];ok=s[3]<2<=s[4];Print[{s[3],s[4],ok}];If[!TrueQ[ok],Exit[1]];
