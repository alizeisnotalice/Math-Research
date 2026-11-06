(* r033: adaptive Bellman computation for the universal PAYOFF, not a CDF.
   State values are linearly interpolated on r>=0.  This is floating point
   evidence only: no interpolation, truncation, or control-error certificate
   is claimed.  The stage kernel is compiled to the Wolfram VM and only the
   previous/current state arrays are retained. *)

ClearAll[compiledStage, valueAt, runCase, makeRows];

mu = N[x /. FindRoot[x - 1 - Exp[-x] == 0, {x, 1.2}], 17];

(* prev[[i]] is V_{k-1}((i-1) h).  The final two returned entries count
   endpoint choices in the c bisection.  They are not errors: r=0,c=mu is
   an exact endpoint root already at the first stage. *)
compiledStage = Compile[{{prev, _Real, 1}, {ntotal, _Integer}, {kstage, _Integer},
    {h, _Real}, {muin, _Real}},
  Module[{m = Length[prev], out, aa, em, aval, i, slope, r, ar, c, clo, chi,
    cm, ee, d0, dhi, dm, x, ax, vx, index, delta, boundary, vneg, aneg,
    endpointLo = 0., endpointHi = 0., iter},
   out = Table[0., {Length[prev] + 2}];
   aval = Table[0., {Length[prev]}];
   em = Exp[-h]; aa = 1. - em;
   boundary = 1. + ((kstage - 1) muin + 1.)/ntotal;
   aval[[1]] = boundary;
   For[i = 1, i < m, i++,
    slope = (prev[[i + 1]] - prev[[i]])/h;
    aval[[i + 1]] = em aval[[i]] + aa prev[[i]] + (h - aa) slope;
    ];
   For[i = 1, i <= m, i++,
    r = (i - 1) h; ar = aval[[i]];
    (* At c=0, A(r-c)=A(r), V(r-c)=V(r). *)
    d0 = 2. (ar - prev[[i]]);
    chi = r + muin;
    ee = Exp[-chi];
    vneg = 1. + kstage muin/ntotal;
    aneg = boundary + muin/ntotal;
    dhi = ar + ee aneg - (1. + ee) vneg;
    If[d0 <= 0.,
     c = 0.; endpointLo = endpointLo + 1.,
     If[dhi >= 0.,
      c = chi; endpointHi = endpointHi + 1.,
      clo = 0.;
      For[iter = 1, iter <= 48, iter++,
       cm = .5 (clo + chi); x = r - cm;
       If[x < 0.,
        ax = boundary - x/ntotal;
        vx = 1. + (kstage - 1) muin/ntotal - x/ntotal,
        index = Floor[x/h] + 1;
        If[index >= m, index = m - 1];
        delta = x - (index - 1) h;
        slope = (prev[[index + 1]] - prev[[index]])/h;
        ee = Exp[-delta];
        ax = ee aval[[index]] + (1. - ee) prev[[index]] +
          (delta - (1. - ee)) slope;
        vx = prev[[index]] + delta slope;
        ];
       ee = Exp[-cm];
       dm = ar + ee ax - (1. + ee) vx;
       If[dm > 0., clo = cm, chi = cm];
       ];
      c = .5 (clo + chi);
      ];
     ];
    x = r - c;
    If[x < 0.,
     ax = boundary - x/ntotal,
     index = Floor[x/h] + 1;
     If[index >= m, index = m - 1];
     delta = x - (index - 1) h;
     slope = (prev[[index + 1]] - prev[[index]])/h;
     ee = Exp[-delta];
     ax = ee aval[[index]] + (1. - ee) prev[[index]] +
       (delta - (1. - ee)) slope;
     ];
    ee = Exp[-c];
    out[[i]] = (ar + ee ax)/(1. + ee);
    ];
   out[[m + 1]] = endpointLo; out[[m + 2]] = endpointHi; out
   ], RuntimeOptions -> "Speed", CompilationTarget -> "WVM"];

valueAt[values_List, x_?NumericQ, h_?NumericQ] := Module[{m = Length[values], j, z},
 If[x <= 0., Return[values[[1]]]];
 j = Min[Floor[x/h] + 1, m - 1]; z = x - (j - 1) h;
 values[[j]] + z (values[[j + 1]] - values[[j]])/h
 ];

runCase[ntotal_Integer, h_?NumericQ] := Module[
  {rmax, m, grid, previous, raw, current, k, t0, qlo, qhi, qm, root,
   monotoneWorst = 0., bracketLow = 0., bracketHigh = 0., stageBadLo = 0.,
   stageBadHi = 0., row},
 rmax = Ceiling[(ntotal mu + 16. Sqrt[ntotal] + 10.)/h] h;
 m = Round[rmax/h] + 1; grid = N[Range[0, m - 1] h];
 (* Correct payoff terminal: V_0(r)=e^-r on r>=0.  The r<0 affine
    continuation is inserted only inside compiledStage. *)
 previous = Developer`ToPackedArray[N[Exp[-grid]]];
 t0 = AbsoluteTime[];
 For[k = 1, k <= ntotal, k++,
  raw = compiledStage[previous, ntotal, k, N[h], mu];
  current = Developer`ToPackedArray[raw[[1 ;; m]]];
  stageBadLo += Round[raw[[m + 1]]]; stageBadHi += Round[raw[[m + 2]]];
  monotoneWorst = Max[monotoneWorst, Max[Differences[current]]];
  previous = current;
  ];
 qlo = 0.; qhi = rmax;
 bracketLow = valueAt[previous, qlo, h] - 1.;
 bracketHigh = valueAt[previous, qhi, h] - 1.;
 If[!(bracketLow >= 0. && bracketHigh <= 0.),
  root = Indeterminate,
  Do[qm = .5 (qlo + qhi);
    If[valueAt[previous, qm, h] >= 1., qlo = qm, qhi = qm], {70}];
  root = .5 (qlo + qhi)
  ];
 row = <|"N" -> ntotal, "h" -> h, "grid_rmax" -> rmax, "grid_points" -> m,
   "q_root_VN_eq_1" -> root, "U_N_at_zero" -> previous[[1]], "V_rmax" -> previous[[-1]],
   "positive_monotonicity_violation" -> monotoneWorst,
   "ideal_arithmetic_control_width" -> (rmax + mu)/2.^48,
   "ideal_arithmetic_q_width" -> rmax/2.^70,
   "root_bracket_low_minus_1" -> bracketLow,
   "root_bracket_high_minus_1" -> bracketHigh,
   "endpoint_c_zero_count" -> stageBadLo, "endpoint_c_upper_count" -> stageBadHi,
   "seconds" -> N[AbsoluteTime[] - t0, 8]|>;
 Print["N=", ntotal, " h=", h, " q=", root, " U_N(0)=", previous[[1]],
  " maxForwardDiff=", monotoneWorst, " endpointCounts=", stageBadLo, "/", stageBadHi,
  " sec=", row["seconds"]]; row
 ];

makeRows[specs_List] := runCase[#[[1]], #[[2]]] & /@ specs;

args = Rest[$ScriptCommandLine];
specs = Which[
   MemberQ[args, "all"],
   Flatten[Table[{nn, hh}, {nn, {16, 32, 64, 128}}, {hh, {.1, .05, .025}}], 1],
   MemberQ[args, "upto64"],
   Flatten[Table[{nn, hh}, {nn, {16, 32, 64}}, {hh, {.1, .05, .025}}], 1],
   MemberQ[args, "n32"], {{32, .1}, {32, .05}, {32, .025}},
   MemberQ[args, "n64"], {{64, .1}, {64, .05}, {64, .025}},
   MemberQ[args, "n128"], {{128, .1}, {128, .05}, {128, .025}},
   True, {{16, .1}, {16, .05}, {16, .025}}];
rows = makeRows[specs];
result = <|"wolfram_version" -> $Version, "method" ->
   "compiled WVM payoff Bellman; exact exponential convolution of piecewise-linear state; continuous c bisection",
  "status" -> "floating point state-interpolation probe only; not certified",
  "terminal" -> "V0(r)=Exp[-r] for r>=0; affine continuation only for r<0",
  "mu" -> mu, "rows" -> rows|>;
out = FileNameJoin[{$TemporaryDirectory, "payoff_bellman33.json"}];
Export[out, result, "RawJSON"];
Print["Wrote ", out];
