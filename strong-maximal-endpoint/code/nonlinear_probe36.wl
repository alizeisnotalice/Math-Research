(* r036: stable machine finite mixtures, checked against r035's 50-digit
   implementation; nonlinear potentials are numerical probes, not certificates. *)
here36=DirectoryName[$InputFileName];args36=$ScriptCommandLine;
out36=If[Length[args36]>=2,args36[[2]],FileNameJoin[{$TemporaryDirectory,"nonlinear_probe36.json"}]];
mode36=If[Length[args36]>=3,args36[[3]],"standard"];
$patch35LibraryMode=True;Get[FileNameJoin[{here36,"patch_probe35.wl"}]];$patch35LibraryMode=False;
mu36=N[mm];p36=1/(1+Exp[mu36]);
Clear[setup36,fast36,case36];
setup36[tt_Integer]:=setup36[tt]=Module[{ks=Range[0,tt],lb,ww},
 lb=N[Log[Binomial[tt,#]]&/@ks,MachinePrecision];ww=Exp[lb-mu36 ks-tt Log[1+Exp[-mu36]]];
 {ks,lb,ww/Total[ww],N[Log[Gamma[tt]],MachinePrecision]}];
fast36[nn_Integer,tt_Integer,rr_?NumericQ]:=fast36[nn,tt,rr]=Module[{ks,lb,ww,lg,x,ii,xx,kk,ll,lmax,pw,ff,fp,fc,fcp,rs,qq,call,uu,up,bc,bp,z,q0},
 If[rr<=0,Return[{1-rr/nn+tt mu36/nn,0.,-1./nn,1-rr/nn+tt mu36/nn,-1./nn,0.,0.,0.}]];
 {ks,lb,ww,lg}=setup36[tt];x=Clip[rr-mu36 ks,{0.,Infinity}];ii=Flatten@Position[(#>0.&)/@x,True];xx=x[[ii]];kk=ks[[ii]];
 ll=lb[[ii]]+(tt-1)Log[xx]-rr-tt Log[1+Exp[-mu36]]-lg;lmax=Max[ll];pw=Quiet[Exp[ll-lmax],General::munfl];ff=Quiet[Exp[lmax+Log[Total[pw]]],General::munfl];pw=pw/Total[pw];
 fp=ff(-1+(tt-1)Dot[pw,1/xx]);fc=ff Dot[pw,xx^2/(tt+1)+mu36^2 kk(tt-kk)/(tt-1)];
 fcp=ff Dot[pw,xx-xx^2/(tt+1)+mu36^2 kk(tt-kk)(1/xx-1/(tt-1))];rs=ff Dot[pw,xx]/tt;
 qq=Quiet[Dot[ww,GammaRegularized[tt,x]],General::munfl];call=Quiet[Dot[ww,tt GammaRegularized[tt+1,x]-(rr-mu36 ks)GammaRegularized[tt,x]],General::munfl];
 uu=p36/8(ff(rr-tt mu36)^2/tt+(N[HarmonicNumber[tt]]-1)fc/tt);
 up=p36/8((fp(rr-tt mu36)^2+2ff(rr-tt mu36))/tt+(N[HarmonicNumber[tt]]-1)fcp/tt);
 (* rho_N differs from N by O(N e^-N); exact 50-digit rho retained here. *)
 z=N[rr-(nn+ProductLog[-nn Exp[-nn]])];
 bc=Quiet[If[z<=0,0.,Exp[-rr+tt Log[z]-LogGamma[tt+1]]+(rr/nn-1)GammaRegularized[tt,0.,z]-tt/nn GammaRegularized[tt+1,0.,z]],General::munfl];
 bp=Quiet[If[z<=0,0.,GammaRegularized[tt,0.,z]/nn-Exp[-rr+tt Log[z]-LogGamma[tt+1]]],General::munfl];
 {rs+qq+call/nn,uu,-rs-qq/nn,1-rr/nn+tt mu36/nn+bc,-1./nn+bp,ff,up,fp}];
validation36=Map[Function[zz,Module[{nn=zz[[1]],tt=zz[[2]],rr=zz[[3]],aa,bb},aa=stats[nn,tt,rr];bb=fast36[nn,tt,rr];
 <|"N"->nn,"t"->tt,"r"->rr,"maximum_absolute_difference"->Max[Abs[N[aa]-bb]],"density_relative_difference"->Abs[N[aa[[6]]]/bb[[6]]-1]|>]],
 {{256,16,18.},{256,16,38.275},{256,64,80.35},{256,64,113.525},{256,128,175.},{256,256,300.}}];
(* Piecewise linear action kernel is inherited from r033. *)
weight36[tt_Integer]:=weight36[tt]=Total[N[Log[Range[2,tt]]^4/Range[2,tt]^(3/2)]];
case36[nn_,tt_,alpha_,beta_,hh_,kind_]:=Module[{grid,st,pr,q,qp,iv,pv,der,active,cv,ans,res,idx,viol,core,et,ep},
 grid=N[Range[0,Ceiling[(tt mu36+5Sqrt[tt mu36])/hh]]hh];st=fast36[nn,tt,#]&/@grid;pr=fast36[nn,tt-1,#]&/@grid;
 et=If[kind==="exponential_summable",weight36[tt],N[HarmonicNumber[tt,2]]];ep=If[kind==="exponential_summable",weight36[tt-1],N[HarmonicNumber[tt-1,2]]];
 q=alpha st[[All,2]]+beta et st[[All,6]];qp=alpha pr[[All,2]]+beta ep pr[[All,6]];
 If[kind==="linear",iv=st[[All,1]]+q;pv=pr[[All,1]]+qp;der=st[[All,3]]+alpha st[[All,7]]+beta et st[[All,8]],
 iv=st[[All,1]]Exp[q/st[[All,1]]];pv=pr[[All,1]]Exp[qp/pr[[All,1]]];
 der=Exp[q/st[[All,1]]]((1-q/st[[All,1]])st[[All,3]]+alpha st[[All,7]]+beta et st[[All,8]])];
 active=Flatten@Position[MapThread[Less,{iv,st[[All,4]]}],True];cv=MapThread[Min,{iv,st[[All,4]]}];pv=MapThread[Min,{pv,pr[[All,4]]}];
 viol=Max[0.,Max[Differences[cv]],Max[Differences[pv]]];ans=compiledStage[pv,nn,tt,hh,mu36];res=Take[ans,Length[cv]]-cv;idx=First@Ordering[res,-1];
 core=Flatten@Position[(#<=2Sqrt[tt Log[tt]]&)/@Abs[grid-tt mu36],True];
 <|"N"->nn,"t"->tt,"kind"->kind,"alpha"->alpha,"beta"->beta,"h"->hh,"maximum_residual"->res[[idx]],"at_r"->grid[[idx]],"residual_over_density"->If[st[[idx,6]]>0,res[[idx]]/st[[idx,6]],Null],
 "core_maximum_residual"->Max[res[[core]]],"grid_positive_increment"->viol,"max_active_inner_derivative"->Max[der[[active]]],
 "interpretation"->"left-screen clip only; PWL state interpolation; not certified"|>];
rows36=Switch[mode36,
 "summable",Flatten[Table[case36[256,tt,1.1,bet,.025,"exponential_summable"],{tt,{16,64,128}},{bet,{.01,.05}}]],
 "refine",{case36[256,128,1.1,.01,.0125,"exponential_summable"]},
 "terminal",{case36[256,256,1.1,.01,.025,"exponential_summable"]},
 _,Flatten[Table[case36[256,tt,1.1,.1,hh,kind],{tt,{16,64,128}},{hh,{.05,.025}},{kind,{"linear","exponential"}}]]];
Export[out36,<|"wolfram"->$Version,"status"->"stable-machine diagnostic, 50-digit spot checks; tiny probabilities may underflow to zero; no interval certificate","validation"->validation36,"rows"->rows36|>,"RawJSON"];
Print[<|"output"->out36,"cases"->Length[rows36],"maximum_spot_difference"->Max[validation36[[All,"maximum_absolute_difference"]]]|>];
