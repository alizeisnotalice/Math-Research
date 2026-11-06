(* Exact checks; the pointwise and geometric claims have analytic proofs. *)
Clear[l,t];checks=<||>;
closed[n_,l_]:=Exp[l]Sum[(-1)^(n-1-j)l^j/j!,{j,0,n-1}]-(-1)^(n-1);
Do[checks["cube_integral_derivative_n"<>ToString[n]]=TrueQ[FullSimplify[D[closed[n,l],l]-Exp[l]l^(n-1)/Factorial[n-1]]===0];checks["cube_integral_origin_n"<>ToString[n]]=TrueQ[(closed[n,l]/.l->0)===0],{n,1,8}];
gm[t_]:=Piecewise[{{0,t<=1},{(t-1)/t,1<t<=2},{1/t,t>2}}];
gp[t_]:=Piecewise[{{0,t<=2},{(t-2)/t,2<t<=4},{2/t,t>4}}];
vals={MaxValue[{gm[t],t>0},t],MaxValue[{gp[t],t>0},t],MaxValue[{(gm[t]+gp[t])/2,t>0},t]};
checks["slab_exact_values"]=TrueQ[vals=={1/2,1/2,3/8}];
If[!And@@Values[checks],Print[checks];Exit[2]];
res=<|"wolfram"->$Version,"status"->"exact symbolic corroboration of slab calibration and finite cube formulas only","checks"->checks,"slab_values"->ToString[vals,InputForm],"all_pass"->True|>;
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/shared_checks55.json"];Export[out,res,"RawJSON"];Print[res];
