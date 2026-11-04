(* Exact identities used by the analytic three-scale proof. *)
Clear[nn, tt, ss, cc, rr, eta, xx];
Print["Wolfram version: ", $Version];
assum = nn > 0 && 0 < tt < nn && ss > 0 && cc > 0;
hardy = Integrate[xx^(tt-nn-1), {xx, ss/cc, Infinity}, Assumptions -> assum];
root = nn/(nn-tt);
lower[rr_] := (rr-1)^tt/rr^nn;
guards = {
 TrueQ[FullSimplify[hardy == cc^(nn-tt) ss^(tt-nn)/(nn-tt), assum]],
 TrueQ[FullSimplify[D[Log[lower[rr]], rr] == (nn-(nn-tt)rr)/(rr(rr-1)), assum && rr>1]],
 TrueQ[FullSimplify[lower[root] == tt^tt (nn-tt)^(nn-tt)/nn^nn, assum]],
 TrueQ[FullSimplify[Limit[tt^tt (nn-tt)^(nn-tt)/nn^nn, tt->0, Direction->"FromAbove"] == 1, nn>0]],
 TrueQ[FullSimplify[Limit[eta^tt (1+eta)^(nn-tt), tt->0, Direction->"FromAbove"] == (1+eta)^nn, nn>0 && 0<eta<1]]
};
Print["Hardy integral: ", hardy];Print["exact guards: ", guards];
If[!And@@guards, Exit[1]];Print["THREE-SCALE GUARDS PASS"];
