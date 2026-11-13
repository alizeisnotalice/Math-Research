out=<|"density_derivative_integral"->ToString[Integrate[t(1+t^2)^(-n/2),{t,0,Infinity},Assumptions->n>2],InputForm],
"critical_slope"->ToString[Solve[D[a s-(1+a^2)/2,a]==0,a],InputForm],
"positive_slope_width_constant"->ToString[N[Log[2 Pi]/2,30],InputForm],
"zero_slope_width_constant"->ToString[N[Log[Pi/2]/2,30],InputForm]|>;
Print[out];Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];
