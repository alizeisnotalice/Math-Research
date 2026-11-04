Print["Independent finite-R review; Wolfram version: ",$Version];
Print["Cancellation used in (4.6):"];
a=(1-2t-u)w;
Print[FullSimplify[2w+(1-2a)/(1-2t)-(1+2w u)/(1-2t)]];
Print["Crude bracket in (4.7), using w<1, delta<u<=1/256, t<=1/65536:"];
bracket=2+(1+2/256)/(1-2/65536);
Print[{bracket,bracket<4}];
Print["Exact rational constant in (4.8):"];
factor=(64/63)(16384/16383)(4096/4095);
Print[{factor,factor<2}];
Print["Lower baseline coefficient:"];
Print[{(63/64)*2^2,(63/64)*2^2>2}];
Print["All-dimensional Bernoulli bookkeeping identities:"];
Print[{
  2n/(128n)==1/64,
  FullSimplify[2n/(128n)^2<=1/16384,Assumptions->n>=2],
  FullSimplify[1/(32n)<=1/64,Assumptions->n>=2]
}];
