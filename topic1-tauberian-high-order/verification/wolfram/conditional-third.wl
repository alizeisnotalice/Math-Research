Clear[n, v, x, r, q, a, rho, u, g, z];
Print["Version: ", $Version];
Print["F10 scalar derivative:"];
Print[FullSimplify[D[(1 + v x^(1/n))^n, x],
  Assumptions -> {n >= 2, 1/2 <= x <= 1, 0 < v < 1}]];
Print["Polynomial capacity divided difference identity residual:"];
Print[FullSimplify[(r - q) Sum[r^(n - 1 - k) q^k, {k, 0, n - 1}] -
  (r^n - q^n), Assumptions -> {n >= 2, Element[n, Integers], r > 0, q > 0}]];
Print["Atypical-mass conversion coefficient:"];
Print[Normal[Series[(1 + v)^n - 1, {v, 0, 3}]]];
Print["Nonempty example v(g), through cubic order, strict slack g^3:"];
vv = 2 g/(1 - 2 g) (1 + g^3)^(1/n);
Print[FullSimplify[Normal[Series[vv, {g, 0, 3}]], Assumptions -> n >= 2]];
Print["Baseline optimal g at that threshold minus actual g, first nonzero:"];
Print[FullSimplify[Limit[(vv/(2 (1 + vv)) - g)/g^4, g -> 0],
  Assumptions -> n >= 2]];
Print["Near-saturation normalized bulk ratio through linear g, treating hole density as q:"];
rr = (1 - 2 g)^n/(1 + g^3);
Print[FullSimplify[Normal[Series[rr, {g, 0, 1}]], Assumptions -> n >= 2]];
Print["Two-layer cancellation: second-source error order n+1 at n>=2:"];
Print[Reduce[n >= 2 && n + 1 < 3, n, Integers]];
Print["Finite rational nonempty checks n=2..12 and g=1/64,1/128,1/256,1/512:"];
checks = Flatten[Table[
  Module[{gg = 1/den, nn = dim, rr0, aa, vv0, capm, ee, kk, norm, tail},
    rr0 = 1 - 2 gg;
    aa = (2 gg/rr0)^nn (1 + gg^3);
    vv0 = aa^(1/nn);
    kk = den^3;
    capm = kk^nn (2 gg)^nn;
    ee = (kk + 2)^nn - capm;
    norm = capm/(aa ee);
    tail = 1 - ((kk - 2)/kk)^nn;
    {nn, den,
      TrueQ[(2 gg)^nn < aa rr0^nn],
      TrueQ[norm >= 1 - 2 (nn + 1) gg],
      TrueQ[tail <= 2 nn/kk],
      TrueQ[gg^3/(1 + gg^3) <= gg^3]}
  ], {dim, 2, 12}, {den, {64, 128, 256, 512}}], 1];
Print[And @@ Flatten[checks[[All, 3 ;; 6]]]];
Print["Number of finite examples checked: ", Length[checks]];
