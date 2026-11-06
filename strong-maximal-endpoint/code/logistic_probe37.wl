(* r037: two-sided bounded uplift; numerical probe only. *)
here37=DirectoryName[$InputFileName];args37=$ScriptCommandLine;
out37=If[Length[args37]>=2,args37[[2]],FileNameJoin[{$TemporaryDirectory,"logistic_probe37.json"}]];
mode37=If[Length[args37]>=3,args37[[3]],"standard"];
$nonlinear36LibraryMode=True;Get[FileNameJoin[{here37,"nonlinear_probe36.wl"}]];$nonlinear36LibraryMode=False;
lift37[l_?NumericQ,v_?NumericQ,u_?NumericQ]:=Module[{a=l-v,q,e},If[u<=0,Return[Min[l,v]]];If[a<=0,Return[l]];If[v<=0,Return[0.]];
 q=Log[v/a]+u l/(v a);If[q>=0,l/(1+Exp[-q]),e=Exp[q];l e/(1+e)]];
case37[nn_,tt_,alpha_,beta_,h_,kind_]:=Module[{grid,st,pr,q,qp,iv,pv,res,ans,idx,core,inc},
 grid=N[Range[0,Ceiling[(tt mu36+5Sqrt[tt mu36])/h]]h];st=fast36[nn,tt,#]&/@grid;pr=fast36[nn,tt-1,#]&/@grid;
 q=alpha st[[All,2]]+beta weight36[tt]st[[All,6]];qp=alpha pr[[All,2]]+beta weight36[tt-1]pr[[All,6]];
 If[kind==="logistic",iv=MapThread[lift37,{st[[All,4]],st[[All,1]],q}];pv=MapThread[lift37,{pr[[All,4]],pr[[All,1]],qp}],
 iv=MapThread[Min,{st[[All,4]],st[[All,1]]Exp[q/st[[All,1]]]}];pv=MapThread[Min,{pr[[All,4]],pr[[All,1]]Exp[qp/pr[[All,1]]]}]];
 inc=Max[0.,Max[Differences[iv]],Max[Differences[pv]]];
 If[inc>10^-12,Return[<|"N"->nn,"t"->tt,"kind"->kind,"h"->h,"positive_increment"->inc,"status"->"nonmonotone; action-root kernel not used"|>]];
 ans=compiledStage[pv,nn,tt,h,mu36];res=Take[ans,Length[iv]]-iv;idx=First@Ordering[res,-1];
 core=Flatten@Position[(#<=2Sqrt[tt Log[tt]]&)/@Abs[grid-tt mu36],True];
 <|"N"->nn,"t"->tt,"kind"->kind,"alpha"->alpha,"beta"->beta,"h"->h,"maximum_residual"->res[[idx]],"at_r"->grid[[idx]],
 "core_maximum_residual"->Max[res[[core]]],"positive_increment"->inc,"status"->"PWL/machine diagnostic; not certified"|>];
rows37=If[mode37==="initial",Flatten[Table[case37[100000000,t,1.1,.01,.025,kind],{t,{16,64,128}},{kind,{"logistic","exponential"}}]],
 Flatten[Table[case37[256,t,1.1,bet,h,kind],{t,{16,64,128}},{bet,{0.,.01}},{h,{.05,.025}},{kind,{"logistic","exponential"}}]]];
Export[out37,<|"wolfram"->$Version,"validation"->validation36,"rows"->rows37|>,"RawJSON"];
Print[<|"output"->out37,"cases"->Length[rows37]|>];
