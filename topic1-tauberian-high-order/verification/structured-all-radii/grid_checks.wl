(* Exact R018 identities in Wolfram Language, all symbols local to this script. *)
ClearAll["Global`*"];
Print["Version ", $Version];
check[label_, truth_] := If[TrueQ[truth], Print[label, " = True"], Print["FAILED ", label, " = ", truth]; Exit[1]];

(* Two inward branches at an actual remote witness constrained by the axial neighbour. *)
branch[pp_, qq_] := Integrate[pp eps^3 (qq eps^3 - y)/(4 nn eps^2), {y, 0, qq eps^3}];
check["two-ended wedge coefficient", Factor[2 (branch[pp, qq]+branch[qq, pp])/(4 eps^2) - pp qq (pp+qq) eps^5/(16 nn)] == 0];
check["uniform wedge erosion", FullSimplify[Min[a0,b0] <= b0, Assumptions -> a0>=0 && b0>=0]];
check["paired face charge", Reduce[Max[0,za-la+lb-gm,zb-lb+la-gm] < (za+zb)/2-gm, {za,zb,la,lb,gm}, Reals] === False];
check["uncaptured mass sign", Reduce[0 < eps < 1/8 && 0 <= ms <= 1/4 && 1-ms-2 eps <= 0, {ms,eps}, Reals] === False];
check["fixed positive mass offset relabelling", Expand[(1-2 (w+zo eps^6/2)+zo eps^6)-(1-2w)] == 0];
check["pointwise gain", Reduce[x>=0 && y>=0 && x y (1-(x+y)/8)>64/27, {x,y}, Reals] === False];

(* The signed planar cell area also follows using the cross product of its diagonals. *)
v0 = {0,0}; v1 = {1+gi eps^6,qq eps^3}; v2 = {pp eps^3,1+gj eps^6};
v3 = {1+hi eps^6+pp eps^3,1+hj eps^6+qq eps^3};
uv = (v1-v0+v3-v2)/2; vv=(v2-v0+v3-v1)/2;
check["quadrilateral opposite-edge identity", Expand[Det[{v3-v0,v2-v1}]/2-Det[{uv,vv}]] == 0];
(* Exact corner and determinant checks use no floating point. *)
Print["OK"];
