Print["Wolfram version: ", $Version];
Print["All computations below are exact symbolic identities, not geometry certificates."];
Print["Isolated two-dimensional corner contact:"];
corner[z_] := (m+z^2)/(4(r0+z)^2);
Print[Factor[D[corner[z],z]]];
Print[FullSimplify[{D[corner[z],z]/.z->m/r0, corner[m/r0]},Assumptions->{m>0,r0>0}]];
Print["All-dimensional first-contact derivative certificate:"];
bracket=1-d-s-2n d-(n-1)z;
certificate=(1-4n d)+n(d-s)+(n-1)(d+s-z);
Print[Expand[bracket-certificate]];
Print["Far-radius common denominator certificate:"];
far=2-(1+2d)/(1+e)-1/(1-2d);
Print[Factor[Together[far]]];
Print[FullSimplify[Together[far]*(1+e)*(1-2d)-(e-4d(1+e-d))]];
Print["Near-radius numerical rational lower bound, n>=2:"];
Print[3*(31/32)*(3/4)];
Print["Partial filling derivative sign certificate:"];
partial[d_] := (2(1-k d))^n/(1-(2d)^n);
Print[FullSimplify[D[partial[d],d]/partial[d],Assumptions->{n>=2,k>3,d>0,d<1/k}]];
Print[FullSimplify[
 D[partial[d],d]/partial[d]-
 n(2^n d^(n-1)-k)/((1-k d)(1-2^n d^n)),
 Assumptions->{n>=2,k>3,d>0,d<1/k}]];
Print["Baseline radius identity and series for reference (not the research result):"];
Print[Table[
 {nn, With[{uu=1/1000},
   FullSimplify[(uu/(1+2uu))/(1-2uu/(1+2uu))==uu]],
   Series[1/(1-(2u/(1+2u))^nn),{u,0,nn+2}]},
 {nn,2,6}]];
