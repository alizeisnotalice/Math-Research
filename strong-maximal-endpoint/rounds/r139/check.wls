u[z_]:=InverseCDF[NormalDistribution[],1-z];Print[N[{ -Sqrt[2]u[1/Sqrt[2]], Sqrt[2]u[1/4], Sqrt[2]u[1/5]},18]]
