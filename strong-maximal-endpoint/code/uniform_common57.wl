(* r057: exact algebra audit for the common-direction response normalizer. *)
ClearAll["Global`*"];
out = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], FileNameJoin[{DirectoryName[$InputFileName], "uniform_common57.json"}]];
Clear[r, b00, mm, b0, bp];
normalizerRewrite = FullSimplify[r (mm + b0)^2/b00 - r^2 mm^2 == r b0^2/b00 + 2 r mm b0/b00 + r^2 mm^2 (1/(r b00) - 1), Assumptions -> r > 0 && b00 > 0];
crossRewrite = FullSimplify[r Abs[mm] Sqrt[b00 bp]/b00 == r Abs[mm] Sqrt[bp/b00], Assumptions -> r > 0 && b00 > 0 && bp >= 0];
quadraticIdentity = FullSimplify[(Sqrt[r b00]^2 - 1) == r (b00 - 1/r), Assumptions -> r > 0 && b00 > 0];
res = <|"status" -> "r057 symbolic normalizer/Schur algebra audit", "wolfram" -> $Version, "normalizer_schur_expansion" -> normalizerRewrite, "cross_term_rewrite" -> crossRewrite, "quadratic_normalizer_identity" -> quadraticIdentity, "scope" -> "Algebra only; analytic response and profile bounds are proved in common_window57.md."|>;
Clear[x, eta];
quadBound=TrueQ[Reduce[0<eta<=1/2&&x>=0&&(x-1)^2<=eta^2 x&&((x-1)>3eta||(x-1)<-3eta),{x,eta},Reals]===False];
res["normalizer_quadratic_bound"]=quadBound;
If[!And@@{normalizerRewrite,crossRewrite,quadraticIdentity,quadBound},Print[res];Exit[2]];
Export[out,res,"RawJSON"];Print[res];
