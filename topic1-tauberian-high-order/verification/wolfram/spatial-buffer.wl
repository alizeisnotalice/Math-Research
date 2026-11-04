(* Exact scalar certificates for the graph-buffer improvement. *)
ClearAll["Global`*"];
SetAttributes[record, HoldRest];
failures = 0;
record[label_String, expression_] := Module[{result = expression},
  Print[label, " = ", ToString[result, InputForm]];
  If[result =!= True, failures++]
];
noExample[variables_List, conditions_] := Resolve[Exists[variables, conditions], Reals] === False;
Print["Kernel: ", $Version];

record["log radii ratio vanishes at zero and has derivative < 8 for q=C1 delta<=1/8",
  ((Log[(1+2 q)/(1-2 q)] /. q -> 0) == 0) &&
  noExample[{q}, 0 <= q <= 1/8 && D[Log[(1+2 q)/(1-2 q)], q] >= 8]];
record["Geometric-series majorant of Exp[1/4] is < 2",
  Sum[(1/4)^m, {m, 0, Infinity}] == 4/3 && 4/3 < 2];
record["A k-step path puts the terminal half-box in a 4k root cube",
  noExample[{k, q, r1, r2}, k >= 1 && 0 <= q <= 1/8 && 0 < r1 <= 2 && 0 < r2 <= 2 &&
    ((1/2+q) (r1+r2) > 3 || 3 k+r2/2 > 4 k)]];
record["Capping graph-distance weights at one preserves edge Lipschitz constant",
  noExample[{x, y, eta}, x >= 0 && y >= 0 && eta >= 0 && Abs[x-y] <= eta &&
    Abs[Min[1, x]-Min[1, y]] > eta]];
record["Edge-length interpolation lower bound",
  noExample[{r, w, Rr, theta}, 0 <= Rr <= r && w >= 0 && 0 <= theta <= 1 &&
    (1-theta) r+theta (Rr+w) < Rr+theta w]];

(* rho=d_i+d_j, w=C1 rho, z=c_n rho, g=C3 rho=4(w+z+rho); all lengths divided by s. *)
g = 4 (w+z+rho);
record["Noncritical margin and common critical-axis sign, all constants abstracted",
  noExample[{rho, w, z}, 0 < rho && w >= 0 && z >= 0 && w+z+rho <= 1/16 &&
    ((1/2-g+z) >= (1-rho) (1/2-w) || (1/2-g-z) <= 0)]];
record["Quantile tail exceeds the source slack after the radius ratio >=1/2",
  noExample[{delta}, 0 < delta < 1/2 && 4 delta (1-delta) <= delta]];

(* U_*, L_* are upper/lower endpoints of the two source/target intervals. *)
commonGap = (1-ti) (l0-u0) + ti (l1-u1);
mixedGap = (1-tj) l0+tj l1 - ((1-ti) u0+ti u1);
record["Mixed normal gap = common normal gap plus target endpoint drift",
  Expand[mixedGap - commonGap - (tj-ti) (l1-l0)] === 0];
record["Common interpolation of two equally oriented disjoint projections has nonnegative gap",
  noExample[{g0, g1, theta}, g0 >= 0 && g1 >= 0 && 0 <= theta <= 1 && (1-theta) g0+theta g1 < 0]];
record["Unequal parameter loss bounded by the Lipschitz increment times endpoint drift",
  noExample[{eta, d, dx, dt}, eta >= 0 && d >= 0 && Abs[dx] <= d && Abs[dt] <= eta && Abs[dt dx] > eta d]];
record["Square of actual two-cell defect charged linearly at regular threshold",
  noExample[{delta, di, dj}, delta >= 0 && 0 <= di <= delta && 0 <= dj <= delta &&
    (di+dj)^2 > 2 delta (di+dj)]];

kap = 1/(n+2);
record["All three exponents in B11 equal 1+kappa",
  Simplify[{1+kap, 2-(n+1) kap, 1+kap} == ConstantArray[1+kap, 3], n >= 2]];
record["Budget in original t and absorption of fixed-n baseline/lower errors",
  Simplify[(n+1+kap)/n == 1+(1+kap)/n && 1+kap < 2 <= n, n >= 2]];
k0 = 1/(2 n+3); d0 = (n+2)/(2 n+3);
record["Coarse two-coordinate fallback, including k delta -> 0",
  Simplify[2*d0 == 1+k0 && 2-n k0-d0 == 1+k0 && d0 > k0, n >= 2]];
Print["Failures: ", failures];
Exit[If[failures === 0, 0, 1]];
