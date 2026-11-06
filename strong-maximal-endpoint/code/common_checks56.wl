$MaxExtraPrecision = 1000;
Print["Wolfram version: ", $Version];
Clear[c, mu, p, m, b, e1, e2, s1, s2, k1, k2, m2];
p[c_] := 1/(1 + Exp[c]);
m[c_] := 1 + c p[c];
(* At exp(-mu)=mu-1, p(mu)=(mu-1)/mu.  Express derivatives only through p. *)
pp = (mu - 1)/mu;
mprimeInP = p0 - c0 p0 (1-p0);
msecondInP = -2 p0 (1-p0) + c0 p0 (1-p0)(1-2p0);
check1 = TrueQ[FullSimplify[(mprimeInP /. {c0->mu, p0->pp}) == 0, Assumptions -> mu>1]];
check2 = TrueQ[FullSimplify[(msecondInP /. {c0->mu, p0->pp}) == -pp, Assumptions -> mu>1]];
rPlus = e1 + e2 + 2 s1 s2 + b m2 (2-k1-k2);
rMinus = e1 + e2 - 2 s1 s2 + b m2 (2-k1-k2);
Print[<|"mPrimeMuZero" -> check1, "mSecondMuMinusP" -> check2,
 "plusDecomposition" -> TrueQ[FullSimplify[rPlus == b m2 (2-k1-k2)+(e1+e2+2s1 s2)]],
 "minusDecomposition" -> TrueQ[FullSimplify[rMinus == b m2 (2-k1-k2)+(e1+e2-2s1 s2)]]|>];

Clear[bb0,bpp,b0p,qq,zz,mm,ee];
psdCheck=TrueQ[FullSimplify[0<=bpp-b0p^2/bb0<=bpp,Assumptions->bb0>0&&bpp>=0&&b0p^2<=bb0 bpp]];
negativeCheck=TrueQ[Reduce[qq>0&&zz>=1/2&&mm<=-3qq/4&&ee<=qq/8&&zz mm+ee> -qq/4,{qq,zz,mm,ee},Reals]===False];
checks=<|"mPrimeMuZero"->check1,"mSecondMuMinusP"->check2,"plusDecomposition"->TrueQ[Expand[rPlus-b m2(2-k1-k2)-(e1+e2+2s1 s2)]===0],"minusDecomposition"->TrueQ[Expand[rMinus-b m2(2-k1-k2)-(e1+e2-2s1 s2)]===0],"psd_response_bounds"->psdCheck,"fixed_input_negative_bound"->negativeCheck|>;
Clear[gg,hh,vv,delta];
quartic=gg delta vv^2/2+hh delta^2/2;
checks["angular_IFT_gain_coefficient"]=TrueQ[FullSimplify[(quartic/.delta->-gg vv^2/(2hh))+gg^2 vv^4/(8hh),Assumptions->hh<0]===0];
If[!And@@Values[checks],Print[checks];Exit[2]];
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/common_checks56.json"];
Export[out,<|"wolfram"->$Version,"status"->"exact symbolic corroboration; asymptotic proof is analytic","checks"->checks,"all_pass"->True|>,"RawJSON"];
