ClearAll[av,aa,qq,vv,mm];
Print["Wolfram ",$Version,"; actual process ",$ProcessID];
Print["R49 exact reverse ratio only; av-aa(1+qq)!=0, aa>0, vv>0. Geometry and asymptotic premises are separate."];
rr=av/(av-aa(1+qq));
lhs=(1+qq)(vv+mm)-av mm/aa;
rhs=(av-aa(1+qq))/aa vv (rr-(1+mm/vv));
Print["exact reverse-ratio residual: ",Together[lhs-rhs]];
Print["terminal positive weight uses Omega>=1: no numerical limit or Tauberian premise is checked."];
Print["complete"];
