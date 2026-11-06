(* Exact symbolic check of the central G cancellations, using exp(mu)=1/(mu-1). *)
Clear[c, d, x, mu, p, u, m, v, fc, fd, dc, m2, gg, red];
p[x_] := 1/(1 + Exp[x]);
u[x_] := x p[x];
m[x_] := 1 + u[x];
v[x_] := m[x] - x D[m[x], x];
fc[c_, d_] := (D[p[x], x] /. x -> c) + u[d]/c^2;
fd[c_, d_] := -(D[u[x], x] /. x -> d)/c;
dc[c_, d_] := d m[d]/(c v[d]);
m2[z_] := D[m[x], {x, 2}] /. x -> z;
gg[c_, d_] := fc[c, d] + fd[c, d] dc[c, d] - m2[c]*(1 - d/c);
red[expr_] := FullSimplify[Together[expr /. {Exp[r_. mu] :> (mu - 1)^(-r),
    E^(r_. mu) :> (mu - 1)^(-r)}],
  Assumptions -> mu > 1];
Print["G(c,c)+m'(c)^2/v(c): ",
 TimeConstrained[FullSimplify[Together[gg[c, c] + (D[m[x], x] /. x -> c)^2/v[c]], Assumptions -> c > 0], 12, "timeout"]];
Print["G(mu,mu): ", TimeConstrained[red[gg[mu, mu]], 12, "timeout"]];
Print["G_c(mu,mu): ", TimeConstrained[red[D[gg[c, d], c] /. {c -> mu, d -> mu}], 12, "timeout"]];
Print["G_d(mu,mu): ", TimeConstrained[red[D[gg[c, d], d] /. {c -> mu, d -> mu}], 12, "timeout"]];
muN = x /. FindRoot[x == 1 + Exp[-x], {x, 1.3}, WorkingPrecision -> 50];
Print["raw symbolic G_c evaluated at numeric root: ", N[D[gg[c, d], c] /. {c -> muN, d -> muN}, 25]];
Print["raw symbolic G_d evaluated at numeric root: ", N[D[gg[c, d], d] /. {c -> muN, d -> muN}, 25]];
