(* Reproduce the exact residual printed by the wolframscript -code audit. *)
mu = 1 + r;
phi = (1 + r Exp[I t mu])/(mu (1 - I t));
pp = r/mu;
a = ComplexExpand[Abs[phi]^2];
b = pp ComplexExpand[Abs[phi - Exp[I t mu]]^2];
Print[FullSimplify[
  1 - a - b - (mu t^2 + 2 r t Sin[mu t] + 2 r^2/mu (1-Cos[mu t]))/(mu^2 (1+t^2)),
  Assumptions -> 0 < r < 1/2 && Element[t, Reals]]];
