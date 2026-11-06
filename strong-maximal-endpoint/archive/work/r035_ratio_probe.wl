(* Probe whether U_t/a_t is decreasing in r below rho_N for the exact
   centered quadratic Poisson formula.  This is evidence only, not a proof. *)
$MaxExtraPrecision = 50000;
mu = N[x /. FindRoot[x == 1 + Exp[-x], {x, 1.3}, WorkingPrecision -> 100], 80];
p = N[1/(1 + Exp[mu]), 80];

polySum[tt_, rr_, pow_, wt_ : (1 &)] := Module[{kk = Range[0, Min[tt, Floor[rr/mu]]]},
 If[Length[kk] == 0, 0``70,
  Total[SetPrecision[Binomial[tt, #] wt[#] (rr-mu #)^pow, 75] & /@ kk]]
];
vals[tt_Integer, nn_Integer, rr_?NumericQ] := Module[
 {kk, ww, uu, ff, aa, ap, cn, cnp, cc, hh, uuval, fprimeoverf, uovera,
  cdf, rrterm, put, at},
 kk = Range[0, Min[tt, Floor[rr/mu]]];
 If[Length[kk] == 0, Return[<|"r" -> rr, "ratio" -> Indeterminate|>]];
 uu = rr-mu kk;
 ww = N[Binomial[tt,#] p^# (1-p)^(tt-#)& /@ kk,70];
 aa = Total[Binomial[tt,#] uu[[# + 1]]^(tt-1)& /@kk];
 ff = Exp[-rr] aa/(1+Exp[-mu])^tt/Gamma[tt];
 cn = polySum[tt,rr,tt+1]/(tt+1)+mu^2/(tt-1) polySum[tt,rr,tt-1,Function[k,k (tt-k)]];
 cc = cn/aa; hh = HarmonicNumber[tt]-1;
 uuval = ff p/8 ((rr-tt mu)^2/tt+cc hh/tt);
 cdf = Total[MapThread[#1 GammaRegularized[tt,0,#2]&, {ww,uu}]];
 rrterm = ff polySum[tt,rr,tt]/(tt aa);
 put = Total[MapThread[#1 (#2 GammaRegularized[tt,0,#2]-tt GammaRegularized[tt+1,0,#2])&, {ww,uu}]];
 at = cdf-rrterm-put/nn;
 <|"t"->tt,"N"->nn,"r"->N[rr,14],"a"->N[at,14],"U"->N[uuval,14],"U/a"->N[uuval/at,14]|>
];

run[tt_,nn_] := Module[{rho, grid, out},
 rho = N[nn + ProductLog[0,-nn Exp[-nn]],60];
 grid = N[Join[Subdivide[1/20, Min[rho-1/10,tt mu-4 Sqrt[tt]],30],
     Subdivide[Max[1/20,tt mu-4 Sqrt[tt]], Min[rho-1/10,tt mu+4 Sqrt[tt]],40]],50];
 out = vals[tt,nn,#]& /@ DeleteDuplicates[grid];
 Export["/tmp/r035_ratio_t"<>ToString[tt]<>".json",out,"RawJSON"];
 Print["t=",tt,"  max successive increase U/a = ",
  Max[Differences[Lookup[out,"U/a"]]],"  endpoints=",Lookup[{First[out],Last[out]},"U/a"]];
];
Scan[run[#,256]&, {16,32,64}];
