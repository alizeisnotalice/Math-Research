(* Wolfram Language symbolic cross-check for the critical left expansion. *)
Print["Wolfram version: ", $Version];
Clear[aSum, q, beta, lambdaC, rhoC, u, cn];
rhoC = aSum/(q - 1);
pLeft = lambdaC (1 + beta/(q - 1))^aSum;
Print["left derivative:", FullSimplify[
  ((-D[pLeft, beta] /. beta -> 0) == -rhoC lambdaC),
  Assumptions -> aSum > 0 && q > 1 && lambdaC > 0]];
Print["inverse Laurent/Taylor series:", Series[
  (q - 1) ((1 + u/lambdaC)^(1/aSum) - 1), {u, 0, 3}]];
super = 1 + cn (1 - beta)/beta;
Print["residue limit:", FullSimplify[
  Limit[(pLeft - lambdaC) super, beta -> 0, Direction -> "FromAbove"],
  Assumptions -> aSum > 0 && q > 1 && lambdaC > 0 && cn > 0]];

(* A critical active-set switch: weights {1,1,2}, rho_c=2. *)
pSwitch = ((1 + beta)/2)^2;
Print["switch case P(1-beta):", pSwitch];
Print["switch derivative:", D[pSwitch, beta] /. beta -> 0];
