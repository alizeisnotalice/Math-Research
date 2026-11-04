(* Independent high-precision / symbolic guard for the cube shell formula. *)
$MaxExtraPrecision=2000;
ClearAll[ip, mu, rr, guard, t];
guard[n_,q_]:=If[TrueQ[q],Print["PASS ",n],Print["FAIL ",n,": ",q];Exit[1]];
ip[lo_,hi_,e_]:=(hi^(e+1)-lo^(e+1))/(e+1);
mu[n_,b_,t_]:=Module[{r=1+t,u=t,p0,p1,p2},
 p0=(1+b)2^(n-1)n/t*u^t;
 p1=2^(n-1)(n ip[u,r,t-1]+b u (n-1) ip[u,r,t-2]);
 p2=(2r)^(n-1)ip[r,r+1,t-n];p0+p1+p2];
rr[n_,b_,t_]:=mu[n,b,t]/(2^n(1+t)^n);
Do[guard["limit_n"<>ToString[n]<>"_b"<>ToString[b],
 FullSimplify[Limit[t rr[n,b,t],t->0,Direction->"FromAbove"]==n(1+b)/2]],
 {n,{1,2,3}},{b,{1,2,5}}];
taus={1/10,1/100,1/10000,1/1000000};
rows=Flatten[Table[{n,b,N[tau rr[n,b,tau],50],N[n(1+b)/2,50]},
 {n,{1,2,3}},{b,{1,2,5}},{tau,taus}],2];
guard["all_numeric_positive",And@@(#[[3]]>0&/@rows)];
Export[FileNameJoin[{DirectoryName[$InputFileName],"cube_shell_wolfram_rows.json"}],
 Association["version"->"round27-cube-shell-wl-v1-20261105","rows"->(AssociationThread[{"n","b","tauR","target"},#]&/@rows)],"RawJSON"];
Print["PASS cube shell guards"];
