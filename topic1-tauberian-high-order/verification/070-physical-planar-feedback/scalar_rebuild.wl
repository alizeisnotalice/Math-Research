Print["Version: ", $Version];
Print["ProcessID: ", $ProcessID];
g = 1/(1 + 2 u);
aa = 4 u^2;
ee = kap u^3;
ww = u g + lam u^6;
rr = 1 - 2 ww;
budget = 4 ww^2 - aa rr^2;
Print["NormalizedBudgetResidual: ",
  Factor[budget/(aa u^5) - (2 lam + (1 - 4 u^2) lam^2 u^5)]];
psi = 2 del z/(1 + Sqrt[1 + 4 a del^2 z]);
Print["InverseResidual: ",
  FullSimplify[psi/del + a psi^2 - z,
    Assumptions -> a > 0 && del > 0 && z >= 0]];
Print["InverseSeries: ", Series[psi, {del, 0, 6}]];
areaPerturbation = 2 Integrate[-a del^3 z^2, {z, 0, 1}];
physicalPerturbation = 2 ee^2 areaPerturbation;
leadingDU = Factor[physicalPerturbation /. {a -> aa, del -> ee/(2 aa rh)}];
leadingDLambda = Factor[leadingDU/(2 aa u^5)];
leadingDM = Factor[(-4 aa rh u^6) leadingDLambda];
coeff = Factor[leadingDM/(aa (2 u)^10)];
Print["OneCornerAreaPerturbation: ", areaPerturbation];
Print["PhysicalBadVolumePerturbation: ", leadingDU];
Print["RootPerturbationLeadingSlope: ", leadingDLambda];
Print["MassGainLeadingDerivative: ", leadingDM];
Print["HaloCoefficientBeforeRhoLimit: ", coeff];
Print["HaloCoefficientRhoLimit: ", coeff /. rh -> 1];
Print["Kappa16over3Coefficient: ", coeff /. {rh -> 1, kap -> 16/3}];
Print["LimitCoefficientResidual: ", Factor[(coeff /. rh -> 1) - kap^5/196608]];
Print["NumericRationalCoefficientResidual: ",
  Factor[(coeff /. {rh -> 1, kap -> 16/3}) - 16/729]];
Print["Scope: exact scalar identities only; no geometric or uniform remainder certification."];
