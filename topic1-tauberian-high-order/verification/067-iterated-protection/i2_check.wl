ClearAll[delta];
data=Get[FileNameJoin[{DirectoryName[$InputFileName],"i2_results.wl"}]];I2=data["I2"];
I1=((delta-1)delta(-100+95delta-60delta^2+26delta^3+4delta^4))/60;
Print["Wolfram ",$Version,"; PID ",$ProcessID];
Print["c3 coefficient in I2-I1: ",Coefficient[I2-I1,delta,3]];
Print["I2 exact samples at delta=1/20,1/40,1/80: ",I2/.delta->#& /@ {1/20,1/40,1/80}];
Print["I2 polynomial degree: ",Exponent[I2,delta]];
