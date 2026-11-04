(* High-precision numerical diagnostics for general p; not a proof certificate. *)
Print["wolfram_version=",$Version];
ps={11/10,3/2,19/10,2,21/10,3,4,8,16};wp=260;
g[p_,d_,q_]:=q (1-d/(2q))^p+(1-q) (1+d/(2(1-q)))^p-1;
gp[p_,d_,q_]:=Module[{del=d/2,a=1-d/(2q),b=1+d/(2(1-q))},a^p-b^p+p del (a^(p-1)/q+b^(p-1)/(1-q))];
rows={};ok=True;bad={};
 Do[
 swExact=2/p; dsExact={Min[swExact/8,1/100],swExact/2,swExact(1-10^-12),swExact,Min[19/10,swExact+Min[1/10,(2-swExact)/2]]};sw=N[swExact,wp]; ds=DeleteDuplicates[N[dsExact,wp]];
 Do[
  If[d<sw-10^-200,
   lo=d/2;hi=1-N[10^-220,wp];
   Do[mid=(lo+hi)/2;If[gp[p,d,mid]<0,lo=mid,hi=mid],{1000}];q=(lo+hi)/2;H=g[p,d,q];res=Abs[gp[p,d,q]];
   test=TrueQ[res<10^-100]&&TrueQ[q>d/2]&&TrueQ[q<1];If[!test,AppendTo[bad,{p,d,q,res,q-d/2,1-q}]];ok=ok&&test,
   q=d/2;H=(1-d/2)^(1-p)-1;res=0
   ];
  inv=If[d>=sw-10^-200,2-2/(H+1)^(1/(p-1)),d];
  If[d>=sw-10^-200,ok=ok&&TrueQ[Abs[inv-d]<10^-100]];
  ok=ok&&TrueQ[H>=0];
  AppendTo[rows,{p,d,q,H,res,Abs[inv-d]}];,
 {d,ds}],
{p,ps}];
Print["rows=",Length[rows]];
Print["max_inner_residual=",Max[rows[[All,5]]]];
Print["max_outer_inverse_error=",Max[Pick[rows[[All,6]],UnitStep[rows[[All,2]]-2/rows[[All,1]]],1]]];
Print["bad_inner=",bad];
ok=ok&&TrueQ[Length[rows]==45];
Print["diagnostic_pass=",ok];
If[!TrueQ[ok],Exit[1]];
