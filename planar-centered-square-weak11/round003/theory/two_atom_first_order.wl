Print["Wolfram ",$Version," PID ",$ProcessID];
ass=0<t<=1/4 && r>0 && s>0 && r^2==1+t/4 && s^2==t;
Print["gap violation ",Resolve[Exists[{t,r,s},ass && r+s/2>=3/2],Reals]];
Print["ratio <=1 violation ",Resolve[Exists[{t,r},0<t<=1/4 && r>0 && r^2==1+t/4 && 9+2t-5r<=4+t],Reals]];
Print["series ",Series[(9+2t-5Sqrt[1+t/4])/(4+t),{t,0,2}]];
Print["at t=1/4 rational enclosure (r<33/32): ", (33/32)^2>17/16];
