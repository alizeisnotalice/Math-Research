(* Exact rational checks for the continuous-curve budget barrier. *)
a=8/5;tau=6643897513452868/10^14;
vUpper=1/Sum[a^k/k!,{k,0,16}]+32 Sum[(-a)^k/(k!(2k+1)),{k,0,12}];
delta=7179/100-tau;expLower=Sum[delta^k/k!,{k,0,30}];
checks=<|"rational_volume_upper_below_21"->TrueQ[vUpper<21],"rational_exp_lower_above_210"->TrueQ[expLower>210]|>;
strongLower=Sum[(35833/500-tau)^k/k!,{k,0,30}];
checks["rational_strong_exp_lower"] = TrueQ[strongLower>8750/47];
If[!And@@Values[checks],Exit[2]];
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/geometry61.json"];
Export[out,<|"wolfram"->$Version,"checks"->checks,"tau_exact"->ToString[tau,InputForm],"volume_rational_upper"->ToString[vUpper,InputForm],"volume_upper_numeric"->N[vUpper,50],"exp_difference_rational_lower"->ToString[expLower,InputForm],"conditional_score_upper_numeric"->N[tau+Log[210],50],"strong_score_upper_numeric"->N[tau+Log[8750/47],50],"strong_required_budget"->"141/1250","scope"->"If the separate directed budget verifier proves exp(tau)P>=1/10, these rational Taylor bounds prove the fixed full parabola-cone score <71.79. No endpoint-constant upper bound."|>,"RawJSON"];Print[{checks,N[vUpper,20],N[tau+Log[210],20]}];
