Print["Wolfram version: ", $Version];
(* The asymmetry coefficient S for n=1 and tau=1/2. *)
s[t_,a_,b_] := 2^t (a^(1/(1-t))+b^(1/(1-t)))^(1-t)/(a+b);
checks={FullSimplify[s[1/2,1,2]==Sqrt[10]/3],FullSimplify[s[1/2,3,3]==1],FullSimplify[Sqrt[10]/3>1]};
Print["collision coefficient identities: ",checks];If[!And@@checks,Exit[1]];
(* primitive shift bound integral_{-R}^R | |x-d|^(t-1)-|x|^(t-1)| dx <= 4 d^t/t *)
(* split the line at 0,d/2,d; cancellation in two tails gives the exact full-line integral *)
exact=FullSimplify[2 ((d/2)^t/t+(d/2)^t/t),Assumptions->0<t<1&&d>0];
Print["full-line translation L1 formula: ",exact];
Print["three-dimensional exponents: ", Table[{n,tt,n-tt}, {n,1,3},{tt,{1/4,1/2,3/4}}]];
Print["Sgap=",N[Sqrt[10]/3-1,50]];
positiveShift=FullSimplify[d^t/t+Limit[((R-d)^t-R^t+d^t)/t,R->Infinity,Assumptions->0<t<1&&d>0],Assumptions->0<t<1&&d>0];
Print["positive half-line translation norm: ",positiveShift];
If[!TrueQ[FullSimplify[positiveShift==2d^t/t,Assumptions->0<t<1&&d>0]],Exit[1]];
Print["collision bound at tau=1/2: ",FullSimplify[12d^(1/2)/(1/2),Assumptions->d>0]];
