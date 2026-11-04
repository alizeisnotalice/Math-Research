(* Exact/high-precision guards for the unbounded cube angular profile. *)
$MaxExtraPrecision=2000;
ClearAll[guard, ea, small, prim, rat, n, g, t, r, a];
guard[name_,q_]:=If[TrueQ[q],Print["PASS ",name],Print["FAIL ",name,": ",q];Exit[1]];
ea[n_,g_]:=1/n+(n-1)/(n-1+g)(1/(1-g)-1/n);
Do[guard["angular_"<>ToString[n]<>"_"<>ToString[g],
 FullSimplify[ea[n,g]==(n-g)/(n(1-g)),0<g<1]],{n,{2,3}},{g,{1/4,1/2,3/4,9/10}}];
small[a_,n_,g_,t_,r_]:=Module[{d=n-1,q=t+g-1,z=1-g},
 If[q==0,2^d a^z ((1+d Log[r/a])/z+d/z^2),
 2^d ((1-d/q)a^t/t+d r^q a^(1-g)/(q (1-g)))]];
prim[l_,n_,g_,t_,r_]:=If[l<=r,small[l,n,g,t,r],small[r,n,g,t,r]+(2r)^(n-1)(l^(t-n+1)-r^(t-n+1))/(t-n+1)];
rat[n_,g_,t_]:=Module[{rr=1+t},(prim[2+t,n,g,t,rr]+prim[t,n,g,t,rr])/(2^n rr^n)];
rows=Flatten[Table[{nn,gg,N[tt rat[nn,gg,tt],50],N[(nn-gg)/(1-gg),50]},
 {nn,{2,3}},{gg,{1/4,1/2,3/4,9/10}},{tt,{1/10,1/100,1/10000,1/1000000}}],2];
guard["all_atom_rows_positive",And@@(#[[3]]>0&/@rows)];
(* The q=0 branch occurs exactly at (gamma,tau)=(9/10,1/10). *)
guard["log_branch_finite",TrueQ[N[rat[2,9/10,1/10],60]>0]];
Export[FileNameJoin[{DirectoryName[$InputFileName],"unbounded_cube_wolfram_rows.json"}],
 Association["version"->"round28-cube-wl-v1-20261105","rows"->(AssociationThread[{"n","gamma","tauR","target"},#]&/@rows)],"RawJSON"];
Print["PASS unbounded cube guards"];
