ClearAll["Global`*"];
Print["VERSION=",$Version];
(* Directed gaps, exact normal-radius identity. *)
gi=sep+ellp+dlp-hii-radi;
gp=sep+ellp-hii+drh-radp;
checks={FullSimplify[(gi-gp)-(radp-radi+dlp-drh)]===0};
(* Logarithmic side-length path: Hessian and integral Taylor coefficient. *)
AppendTo[checks,FullSimplify[D[Log[1+kap Exp[z]],{z,2}]-kap Exp[z]/(1+kap Exp[z])^2]===0];
(* Scalar affine-layer telescope, with independent positive reserves. *)
steps=5; w[h_]:=tau^(steps-h); tau=bb vv/cc;
expr=Sum[w[h](cc ff[h]-bb vv ff[h-1]),{h,1,steps}];
AppendTo[checks,FullSimplify[expr-(cc ff[steps]-bb vv tau^(steps-1) ff[0])]===0];
(* Exact coefficient on the source/target defect under the root normalization. *)
AppendTo[checks,FullSimplify[vv (1+vv)^(nn-1)-D[(1+vv (1+eps)^(1/nn))^nn,eps]/.eps->0,
 Assumptions->nn>=2&&vv>0]===0];
Print["CHECKS=",checks];Print["FAILURES=",Count[checks,Except[True]]];
(* This is a scale check, not a proof of compactness. *)
Print["CRITICAL_PHYSICAL_ASPECT_L2_SQUARED=",FullSimplify[vv^2 vv^5/vv^6]];
Print["CRITICAL_FACE_JUMP_over_v4_L1=",FullSimplify[vv^6/vv^4]];
Print["CDF_SHIFT_over_v4_L1=",FullSimplify[vv vv^5/vv^4]];

ClearAll[v,k,n,dcell,base];
dcell=(1+v)^(n-2)(1+v Exp[k v^2])(1+v Exp[-k v^2]);base=(1+v)^n;
Print["ANISOTROPY_EXACT_DIFF=",FullSimplify[dcell-base-2v(1+v)^(n-2)(Cosh[k v^2]-1),Assumptions->v>0&&n>=2]];
Print["ANISOTROPY_RATIO_COEFF=",FullSimplify[Limit[-(dcell-base)/(v^5 (dcell-v^n)(base-v^n)),v->0,Direction->"FromAbove",Assumptions->n>=2&&Element[n,Integers]&&Element[k,Reals]]]];
Print["ANISOTROPY_MIXED_NORMALIZED=",FullSimplify[Limit[-v^2(dcell-base)/(v^6 (dcell-v^n)(base-v^n)),v->0,Direction->"FromAbove",Assumptions->n>=2&&Element[n,Integers]&&Element[k,Reals]]]];
