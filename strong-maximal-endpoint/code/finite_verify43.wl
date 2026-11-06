$finite43LibraryMode=True;Get[FileNameJoin[{DirectoryName[$InputFileName],"finite_left43.wl"}]];$finite43LibraryMode=False;
argsV43=$ScriptCommandLine;outV43=If[Length[argsV43]>=2,argsV43[[2]],"/tmp/finite_verify43.json"];
rowsV43={};rawV43={};Do[configure43[tt^power,1/10];ll=N[(fraction/((11/10)mu43(mu43-1)/8 tt))^(1/3),70];x=case43[tt,ll,11/10,64,1/16,True];AppendTo[rowsV43,x];xx=case43[tt,ll,11/10,64,1/16,False];AppendTo[rawV43,xx];Print[KeyTake[x,{"t","N","source","active","residual_scaled_lambda2","seconds"}]],{pars,{{512,2,1/4},{512,8,1},{8192,8,1/4},{8192,8,1}}},{tt,{pars[[1]]}},{power,{pars[[2]]}},{fraction,{pars[[3]]}}];
std43[64]=.;wp43=80;configure43[512^2,1/10];ll=N[((1/4)/((11/10)mu43(mu43-1)/8 512))^(1/3),100];prec43=case43[512,ll,11/10,64,1/16,True];
(* Independent finite mixture includes exactly the u/N terminal defect. *)
tt=512;setup43[tt,ll,1/16];rr=rr43;ks=Range[0,tt];xx=Clip[rr-a43 ks,{0,Infinity}];bi=N[Binomial[tt,#]&/@ks,110];ww=bi Exp[-a43 ks]/(1+Exp[-a43])^tt;
fExact=Exp[-rr]Dot[bi,xx^(tt-1)]/(Gamma[tt](1+Exp[-a43])^tt);cdf=Dot[ww,GammaRegularized[tt,0,xx]];first=Dot[ww,tt GammaRegularized[tt+1,0,xx]+a43 ks GammaRegularized[tt,0,xx]];
dExact=cdf-Exp[-rr]Dot[bi,xx^tt]/(Gamma[tt+1](1+Exp[-a43])^tt)-(rr cdf-first)/N43;
ind43=<|"t"->tt,"N"->N43,"density_relative_difference"->ToString[(fc43-fExact)/fExact,InputForm],"finite_A_over_f_difference"->ToString[dn43-dExact/fExact,InputForm]|>;
Export[outV43,<|"wolfram"->$Version,"status"->"finite matched quadrature refinement, raw sufficient checks, 80-digit repetition and independent mixture; not certified","refinements"->rowsV43,"raw_checks"->rawV43,"precision80"->prec43,"independent_mixture"->ind43|>,"RawJSON"];Print[ind43];
