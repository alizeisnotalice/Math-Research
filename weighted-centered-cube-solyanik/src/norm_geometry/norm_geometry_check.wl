Print["Wolfram version: ",$Version];
Clear[n,tau];
window = FullSimplify[
 Equivalent[tau/(2 n)<(n-2 tau)/(n (n-tau)),
            tau^2-(n+4) tau+2 n>0],
 Assumptions -> n>0 && 0<tau<n];
Print["Window equivalence (positive quadratic): ",window];
If[window =!= True,Exit[1]];
largeDimension = FullSimplify[(n+4-Sqrt[n^2+16])/2>=1,
 Assumptions -> n>=3];
Print["T_n=1 for n>=3: ",largeDimension];
If[largeDimension =!= True,Exit[1]];
Print["T1,T2: ",FullSimplify[{(5-Sqrt[17])/2,(6-Sqrt[20])/2}]];
constant = FullSimplify[
 ((n-tau)*2*n*2^(n+1-tau))/(n/tau) == tau*(n-tau)*2^(n+2-tau),
 Assumptions -> n>0 && 0<tau<n];
Print["Thin-shell constant: ",constant];
If[constant =!= True,Exit[1]];
Print["PASS: three exact symbolic assertions."];
