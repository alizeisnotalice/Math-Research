(* Exact cost and budget/point-regret tradeoff of the r057 actual step inputs. *)
Clear[n,m,del,rho,jbar,s,z];
phi[n_,h_]:=If[h<=1,h,h(1+Log[h])^(n-1)];
checks=<||>;
Do[gm=Table[1-If[i==j,dd,0]+If[i<mm&&j==i+1,dd,0],{i,mm},{j,mm}];
 actual=Total[Map[phi[nn,#]&,gm,{2}],2];
 formula=mm^2-dd+(mm-1)(1+dd)((1+Log[1+dd])^(nn-1)-1);
 checks["actual_Phi"<>ToString[nn]<>"_M"<>ToString[mm]<>"_delta"<>ToString[dd]]=TrueQ[FullSimplify[actual-formula]===0],{nn,{2,4,8,16}},{mm,{2,3,5}},{dd,{1/2,1/3}}];
checks["derivative_lower_from_log_upper"]=TrueQ[FullSimplify[1/((1+s)(1+z))>=1/(1+s)^2,Assumptions->0<=s<=1/2&&0<=z<=s]];
checks["rational_derivative_lower_4over9"]=TrueQ[FullSimplify[1/(1+s)^2>=4/9,Assumptions->0<=s<=1/2]];
checks["normalized_cost_rewrite"]=TrueQ[FullSimplify[(m^2-del+(m-1)(1+del)((1+Log[1+del])^(n-1)-1))/m^2-(1-del/m^2+(m-1)/m^2 (phi[n,1+del]-(1+del))),Assumptions->del>0]===0];
checks["matrix_coefficient"]=TrueQ[FullSimplify[m^2/(m-1)<=2m,Assumptions->m>=2]];
checks["regret_substitution"]=TrueQ[FullSimplify[(4/9)(n-1)del-((4/9)(n-1)m^2 rho)/.del->m^2 rho]===0];
If[!And@@Values[checks],Print[checks];Exit[2]];
(* Exact symbolic input parameters; high-precision numeric diagnostics only. *)
rows=Flatten[Table[dd=Min[1/2,mm^2 Log[nn]^2/nn];jc=(mm^2-dd+(mm-1)(1+dd)((1+Log[1+dd])^(nn-1)-1))/mm^2;
 <|"N"->nn,"M"->mm,"delta"->ToString[N[dd,40],InputForm],"regret"->ToString[N[dd/mm^2,40],InputForm],"log_normalized_budget"->ToString[N[Log[jc],40],InputForm],"log_normalized_budget_80digits"->ToString[N[Log[jc],80],InputForm],"precision_double_absolute_difference"->ToString[Abs[N[Log[jc],80]-N[Log[jc],40]],InputForm],"precision_double_check"->TrueQ[Abs[N[Log[jc],80]-N[Log[jc],40]]<10^-30],"critical_target"->ToString[N[Log[nn]^2/nn,40],InputForm],"delta_hits_cap"->TrueQ[dd==1/2]|>,{nn,{64,256,1024,4096,16384}},{mm,{2,4,8}}],1];
If[!And@@Lookup[rows,"precision_double_check"],Print["Precision doubling failed"];Exit[3]];
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/budget58.json"];
res=<|"wolfram"->$Version,"checks"->checks,"all_pass"->True,"numeric_scope"->"40-digit diagnostics, not directed certificates; exact inequalities proved analytically","rows"->rows|>;
Export[out,res,"RawJSON"];Print[<|"all_pass"->True,"exact_checks"->Length[checks],"rows"->Length[rows],"output"->out|>];
