Clear[tau,z];
Print["Wolfram version: ",$Version];
Print["Window: ",Reduce[0<tau<1 && tau/2<(1-2 tau)/(1-tau),tau,Reals]];
Print["Shift series: ",Series[((1+z)^tau+(1-z)^tau)/2,{z,0,4}]];
Print["Threshold: ",N[(5-Sqrt[17])/2,30]];
