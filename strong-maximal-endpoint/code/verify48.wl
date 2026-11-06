$endpoint48LibraryMode=True;Get[FileNameJoin[{DirectoryName[$InputFileName],"endpoint48.wl"}]];$endpoint48LibraryMode=False;
(* Independent small-n finite truncated powers for log-density and score. *)
finitePoint48[n_,c_,q_,wp_]:=Module[{ks=Range[0,n],pos,terms,summ,sd,lf,score,k},pos=Clip[q-c ks,{0,Infinity}];terms=N[Binomial[n,#]&/@ks,wp];summ=Dot[terms,pos^(n-1)];sd=(n-1)Dot[terms,pos^(n-2)];lf=-q+Log[summ]-LogGamma[n]-n Log[1+Exp[-c]];score=n p46[c]-(n-1)Dot[terms ks,pos^(n-2)]/summ;k=sd/summ;<|"logf"->lf,"score"->score,"k"->k|>];
small48=Table[Module[{row=case48[256,label],c,q,pt,fi,diffs},c=ToExpression[row["c_internal"]];q=ToExpression[row["q_internal"]];pt=point48[256,c,q,120,1/16,32];fi=finitePoint48[256,c,q,120];diffs=Abs[{Log[pt["f"]]-fi["logf"],pt["score"]-fi["score"],pt["log_R_derivative"]-fi["k"]}];<|"label"->label,"logf_score_k_differences"->(ToString[#,InputForm]&/@diffs),"all_less_than_1e-25"->TrueQ[Max[diffs]<10^-25]|>],{label,{"mu","sqrtL","z"}}];
refs48=Table[Module[{a=case48[n,"z"],b=case48[n,"z",2,1/32,40],diff},diff=Abs[a["score_q"]-b["score_q"]];<|"N_exact"->ToString[n,InputForm],"score_difference"->ToString[diff,InputForm],"less_than_1e-20"->TrueQ[diff<10^-20],"roots_valid"->TrueQ[a["root_valid"]&&b["root_valid"]]|>],{n,{256,2^48,2^1024}}];
(* Exact leading symmetric cancellation for a generic smooth u at mu. *)
c0=mm-mm hh/2;dd=c0(1+hh+cc hh^2);uT[x_]:=uu-pp(x-mm)^2/2+vv(x-mm)^3/6;
cancel48=FullSimplify[Coefficient[Normal[Series[uT[dd]-uT[c0],{hh,0,2}]],hh,2]==0];
res48=<|"wolfram"->$Version,"symmetric_quadratic_cancellation"->cancel48,"independent_finite_kernel"->small48,"precision_mesh_cutoff_refinements"->refs48,"scope"->"exact polynomial identity and independent floating checks, no total interval certificate"|>;
out48=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/verify48.json"];Export[out48,res48,"RawJSON"];Print[res48];
