(* Exact elementary components for the scoped pair/geometry theorem. *)
Clear[x,a,e,v0,v1,d0,d1,w0,w1,z,s];
checks=<||>;
vol=2-Exp[-2a e];
checks["cone_volume_at_zero"]=TrueQ[(vol/.e->0)==1];
checks["right_log_volume_first"]=TrueQ[FullSimplify[D[Log[vol],e]/.e->0]==2a];
checks["right_log_volume_second"]=TrueQ[FullSimplify[D[Log[vol],{e,2}]/.e->0]==-8a^2];
(* Segment width dx=(v1-v0)+(d1-d0); z=dx+dy=2(d1-d0). *)
z=2(d1-d0);f[z_]:=(1-Exp[-z])/z;
segment=(v1-v0+d1-d0)Exp[-2d0]f[z];
angular=(v1-v0)Exp[-2d0]f[z];
checks["segment_angular_boundary_decomposition"]=TrueQ[FullSimplify[segment-angular-(Exp[-2d0]-Exp[-2d1])/2,Assumptions->d1!=d0]===0];
checks["F_integral"]=TrueQ[FullSimplify[Integrate[Exp[-s x],{x,0,1}]-(1-Exp[-s])/s,Assumptions->s!=0]===0];
checks["F_removable_limit"]=TrueQ[Limit[(1-Exp[-s])/s,s->0]===1];
(* Envelope bookkeeping uses exact formal differentiation. *)
Clear[la,tt,bb,pp];
checks["pair_minus_single_response"]=TrueQ[Expand[(-la tt)-(-la(tt+bb))]===la bb];
If[!And@@Values[checks],Print[checks];Exit[2]];
res=<|"wolfram"->$Version,"status"->"exact symbolic corroboration only; geometric/variational hypotheses proved separately","checks"->checks,"all_pass"->True,"log_volume_right_series"->ToString[Series[Log[vol],{e,0,2}],InputForm]|>;
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/pair_algebra55.json"];
Export[out,res,"RawJSON"];Print[res];
