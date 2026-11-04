(* Independent exact baseline and singular-central-packet checks. *)
Print["Wolfram version: ", $Version];
Do[
  haloMass = 2 (5/3)^tau/tau;
  inputCost = 8/3;
  lower = Together[haloMass/inputCost];
  expected = Together[3/(4 tau) (5/3)^tau];
  Print["tau=", tau, ", pair identity=", FullSimplify[lower==expected],
        ", lower=", lower],
  {tau,{1/4,1/2}}
];
Do[
  (* Unit Lebesgue mass on (-eps,eps) has this weighted L1 cost. *)
  cost = Integrate[Abs[x]^(tau-1)/(2 eps),{x,-eps,eps},
                   Assumptions->eps>0];
  Print["tau=",tau,", central-packet cost=",cost,
        ", divergence=",FullSimplify[Limit[cost,eps->0,Direction->"FromAbove"]==Infinity]],
  {tau,{1/4,1/2}}
]
