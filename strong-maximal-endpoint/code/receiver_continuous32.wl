(* r032: split-branch receiver-parameter checks and uniform-parameter control.
   Floating point evidence only.  Run: wolframscript -file receiver_continuous32.wl *)

base = DirectoryName[DirectoryName[$InputFileName]];
seedPath = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]],
   FileNameJoin[{base, "rounds", "r030", "receiver_seed.json"}]];
outPath = If[Length[$ScriptCommandLine] >= 3, $ScriptCommandLine[[3]],
   FileNameJoin[{$TemporaryDirectory, "receiver_continuous32.json"}]];
seed = Import[seedPath, "RawJSON"];
n = seed["n"]; b = n/2; cap = seed["height_cap"]; c = seed["c"]; tau = seed["tau_ref"];
w0 = N[seed["weights"]]; vv = N[seed["v_nonnegative"]]; dd = N[seed["receiver_d"]];
ng = Length[w0]; coeff = Join[{1.}, ConstantArray[2., ng - 1]];
u = Join[{{0., 0.}}, Transpose[{Rest[dd] + Rest[vv], Rest[dd] - Rest[vv]}],
   Transpose[{Rest[dd] - Rest[vv], Rest[dd] + Rest[vv]}]];
p = 1./(1. + Exp[c]); mz = 1. + c p; vz = 1. + c^2 p (1. - p);
lo = Min[u]; hi = Max[u] + b mz + 16. Sqrt[b vz];
baseScale = 1.0486657924183962; target = 1.00001;
ks = Range[0, b]; bin = N[Binomial[b, #] & /@ ks]; polynorm = 1./(Gamma[b] (1. + Exp[-c])^b);
gbw = bin Exp[-c ks]/(1. + Exp[-c])^b;

points[j_Integer] := If[j == 1, {{0., 0.}}, {{u[[j, 1]], u[[j, 2]]}, {u[[j + ng - 1, 1]], u[[j + ng - 1, 2]]}}];

buildPolyCV[bb_Integer,cx_Real]:=With[{cc=cx,degree=bb-1,bin=Developer`ToPackedArray[N[Binomial[bb,#]&/@Range[0,bb]]],norm=N[1/(Gamma[bb](1+Exp[-cx])^bb)],nk=bb+1},Compile[{{tt,_Real,1}},Module[{res=Table[0.,{Length[tt]}],i=0,k=0,x=0.,s=0.},For[i=1,i<=Length[tt],i++,s=0.;For[k=1,k<=nk,k++,x=tt[[i]]-cc(k-1);If[x>0.,s=s+bin[[k]]x^degree]];res[[i]]=norm s];res],RuntimeOptions->"Speed"]];

pc = buildPolyCV[b, c];
polyVec[tt_List] := pc[Developer`ToPackedArray[N[tt]]];
gbcdf[z_?NumericQ] := Total[gbw GammaRegularized[b, 0., Clip[z - c ks, {0., Infinity}]]];
phi[h_?NumericQ] := h (1. + Max[Log[Max[h, 1.]], 0.])^(n - 1);
phiPrime = (1. + Log[cap])^(n - 2) (n + Log[cap]);

(* One receiver-density vector and likelihood vector.  The latter is the
   likelihood relative to exp(-t1-t2) after the tau_ref scaling. *)
extendWeights[ww_] := Join[{ww[[1]]}, Rest[ww], Rest[ww]];
sourceAmp[t1_?NumericQ, ww_] := extendWeights[ww] polyVec[t1 - u[[All, 1]]] Exp[Total[u, {2}] - tau];
likeAmp[t2_?NumericQ, amp_] := Dot[amp, polyVec[t2 - u[[All, 2]]]];
profile[s_?NumericQ] := Module[{z, r, d}, If[s < 1., Return[0.]]; If[s <= n, Return[1.]];
  d = Log[s]; z = Max[Exp[(d - Log[n])/(n - 1)] - 1., 0.];
  Do[r = (n - 2) Log[1. + z] + Log[n + z] - d; z = Max[0., z - r/((n - 2)/(1. + z) + 1./(n + z))], {10}]; Min[Exp[z], cap]];

quadCache = <||>;
stdQuad[ord_Integer] := If[KeyExistsQ[quadCache, ord], quadCache[ord], quadCache[ord] = Module[{j, ev, vec, ix, x, wt},
  j = N[DiagonalMatrix[Range[1, ord - 1]/Sqrt[4 Range[1, ord - 1]^2 - 1], 1] +
      DiagonalMatrix[Range[1, ord - 1]/Sqrt[4 Range[1, ord - 1]^2 - 1], -1]];
  {ev, vec} = Eigensystem[j]; ix = Ordering[ev]; x = ev[[ix]]; wt = 2. vec[[ix, 1]]^2;
  {x, wt}]];
quad[ord_Integer, aa_?NumericQ, bb_?NumericQ] := Module[{q = stdQuad[ord]}, {((bb - aa)/2.) q[[1]] + (aa + bb)/2., ((bb - aa)/2.) q[[2]]}];
outerQuad[nOuter_Integer] := Module[{pan = Ceiling[nOuter/16.], edge, x, wt, outx = {}, outw = {}},
  edge = Subdivide[lo, hi, pan]; Do[{x, wt} = quad[16, edge[[i]], edge[[i + 1]]]; outx = Join[outx, x]; outw = Join[outw, wt], {i, pan}]; {outx, outw}];

boundary[lev_?NumericQ, amp_] := Module[{l = lo, r = hi, mid, fl, fr},
  fl = likeAmp[l, amp] - lev; fr = likeAmp[r, amp] - lev;
  If[fl >= 0., Return[l]]; If[fr <= 0., Return[r]];
  Do[mid = (l + r)/2.; If[likeAmp[mid, amp] < lev, l = mid, r = mid], {45}]; (l + r)/2.];

(* Receiver-dependent c_j check.  The source likelihood and H remain fixed
   at c; only the probability law used for each receiver average uses cj.
   Thus an observed cost reduction after rescaling is not a fixed-c dual move. *)
polyC[t_?NumericQ, cj_?NumericQ] := If[t <= 0., 0., Total[bin Clip[t - cj ks, {0., Infinity}]^(b - 1)]/(Gamma[b] (1. + Exp[-cj])^b)];
gbcdfC[z_?NumericQ, cj_?NumericQ] := Module[{ww = bin Exp[-cj ks]/(1. + Exp[-cj])^b}, Total[ww GammaRegularized[b, 0., Clip[z - cj ks, {0., Infinity}]]]];
probVecC[t1_?NumericQ, t2_?NumericQ, cc_List] := Module[{ans = ConstantArray[0., ng], pt, z, q1, q2},
 Do[pt = points[j]; Do[q1 = polyC[t1 - z[[1]], cc[[j]]]; q2 = polyC[t2 - z[[2]], cc[[j]]];
   ans[[j]] += q1 q2 Exp[z[[1]] + z[[2]] - t1 - t2], {z, pt}], {j, ng}]; ans];
constantPartC[t1_?NumericQ, left_?NumericQ, right_?NumericQ, h_?NumericQ, cc_List] := Module[{cost = 0., av = ConstantArray[0., ng], pt, z, fac},
 If[right <= left, Return[{cost, av}]]; cost = Exp[tau - t1] (Exp[-left] - Exp[-right]) phi[h];
 Do[pt = points[j]; Do[fac = polyC[t1 - z[[1]], cc[[j]]] Exp[z[[1]] - t1];
   If[fac > 0., av[[j]] += fac (gbcdfC[right - z[[2]], cc[[j]]] - gbcdfC[left - z[[2]], cc[[j]]])], {z, pt}], {j, ng}]; {cost, av}];
softPartC[t1_?NumericQ, left_?NumericQ, right_?NumericQ, ord_Integer, amp_, cc_List] := Module[{x, wt, cost = 0., av = ConstantArray[0., ng], h},
 If[right <= left, Return[{cost, av}]]; {x, wt} = quad[ord, left, right];
 Do[h = profile[likeAmp[x[[i]], amp]]; cost += wt[[i]] Exp[tau - t1 - x[[i]]] phi[h];
   av += wt[[i]] h probVecC[t1, x[[i]], cc], {i, Length[x]}]; {cost, av}];
evaluateReceiverC[cc_List, srcScale_?NumericQ, nOuter_Integer, inner_Integer] := Module[{ww = srcScale w0, tx, tw, cost = 0., av = ConstantArray[0., ng], roots, q1, a, amp, start = AbsoluteTime[]},
 {tx, tw} = outerQuad[nOuter];
 Do[amp = sourceAmp[tx[[i]], ww]; roots = boundary[#, amp] & /@ {1., N[n], phiPrime};
   {q1, a} = constantPartC[tx[[i]], roots[[1]], roots[[2]], 1., cc]; cost += tw[[i]] q1; av += tw[[i]] a;
   {q1, a} = constantPartC[tx[[i]], roots[[3]], hi, cap, cc]; cost += tw[[i]] q1; av += tw[[i]] cap a;
   {q1, a} = softPartC[tx[[i]], roots[[2]], roots[[3]], inner, amp, cc]; cost += tw[[i]] q1; av += tw[[i]] a, {i, Length[tx]}];
 <|"source_scale" -> srcScale, "receiver_c" -> cc, "outer_nodes" -> Length[tx], "soft_inner_order" -> inner,
  "scaled_cost" -> cost, "actual_logJ" -> Log[cost] - tau, "averages_by_group" -> av/coeff,
  "min_average" -> Min[av/coeff], "seconds" -> AbsoluteTime[] - start|>];


(* Changes in c are accepted only after recomputing averages; no FOSD in c. *)
localC=ConstantArray[c,ng];localC[[{23,24,25}]]+=.05;
allC=ConstantArray[c+.05,ng];trialScale=1.0471293307430003;
localAtBase=evaluateReceiverC[localC,baseScale,256,32];
allAtBase=evaluateReceiverC[allC,baseScale,256,32];
local32=evaluateReceiverC[localC,trialScale,256,32];
local64=evaluateReceiverC[localC,trialScale,256,64];
result=<|"wolfram"->$Version,"status"->"floating continuous split-branch diagnostics; no interval certificate; no parameterwise stochastic order","source_c"->c,"base_source_scale"->baseScale,"trial_source_scale"->trialScale,"changed_groups"->{23,24,25},"receiver_c_increment"->.05,"finite_domain"->{lo,hi},"local_at_base"->localAtBase,"all_at_base"->allAtBase,"local_candidate_soft32"->local32,"local_candidate_soft64"->local64|>;
Export[outPath,result,"RawJSON"];Print[<|"local_base_min"->localAtBase["min_average"],"all_base_min"->allAtBase["min_average"],"candidate_scaled_cost"->local64["scaled_cost"],"candidate_min"->local64["min_average"]|>];
