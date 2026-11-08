a82=Rest[$ScriptCommandLine];input82=a82[[1]];out82=a82[[2]];If[!DirectoryQ[out82],CreateDirectory[out82,CreateIntermediateDirectories->True]];
ClearAll[mean82,root82,run82];
mean82[n_Integer,aa_Integer,lam_,q_?NumberQ,wp_Integer]:=Module[{rt=N[Sqrt[n],wp],top,gg,first,ii},
 top=Min[(aa-1)rt,(n-q)/rt+40];gg[w_?NumericQ]:=With[{l=Log[1+w/rt]},q+l+(n-1)Log[1+l]];
 first=Exp[n Log[q]-q-LogGamma[n+1]];
 ii=NIntegrate[GammaRegularized[n,gg[w]]/rt,Evaluate[Join[{w,0},Select[{1,2,4,8,16,32},#<top&],{top}]],WorkingPrecision->wp,AccuracyGoal->Floor[wp/2]+2,PrecisionGoal->Floor[wp/2]+2,MaxRecursion->18,Method->{"GlobalAdaptive","SymbolicProcessing"->0}];
 first-GammaRegularized[n,0,q]+ii+1-lam];
run82[z_]:=Catch[Module[{n=z["N"],aa=z["A"],guess=Rationalize[z["q"],0],rows={},q,ze,p,rt,zg,r,u,tail,cut,lam=1,relax},
 Do[q=x/.FindRoot[mean82[n,aa,lam,x,wp]==0,{x,N[guess-1/10,wp],N[guess+1/10,wp]},WorkingPrecision->wp,AccuracyGoal->Floor[wp/2]-6,PrecisionGoal->Floor[wp/2]-6,MaxIterations->35];
 rt=N[Sqrt[n],wp];zg=x/.FindRoot[Erfc[x/Sqrt[2]]/2==x/rt,{x,N[z["gaussian_z"],wp]},WorkingPrecision->wp,AccuracyGoal->Floor[wp/2],PrecisionGoal->Floor[wp/2]];r=rt Exp[-zg^2/2]/Sqrt[2Pi];u=(r(zg^2-1)/3+r+zg^2)/(r+1);ze=Log[1+Log[aa]]/Log[aa];p=1+ze(n-1);cut=Min[(aa-1)rt,(n-q)/rt+40];tail=If[cut<(aa-1)rt,(aa-1)GammaRegularized[n,q+Log[1+cut/rt]+(n-1)Log[1+Log[1+cut/rt]]],0];
 relax=If[aa==2, x/.FindRoot[mean82[n,aa,1-1/n^2,x,wp]==0,{x,q,q+N[1/n,wp]},WorkingPrecision->wp,AccuracyGoal->Floor[wp/2]-6,PrecisionGoal->Floor[wp/2]-6,MaxIterations->35],q];
 AppendTo[rows,<|"precision"->wp,"q"->ToString[N[q,wp-10],InputForm],"mean_residual"->ToString[N[mean82[n,aa,lam,q,wp],30],InputForm],"integration_tail_upper"->ToString[N[tail,30],InputForm],"log_charge_saving"->ToString[N[n p Log[p/(p-1)]-q,wp-10],InputForm],"gaussian_reference"->ToString[N[n-rt zg,wp-10],InputForm],"asymptotic_root_prediction"->ToString[N[n-rt zg+u,wp-10],InputForm],"root_minus_asymptotic"->ToString[N[q-(n-rt zg+u),40],InputForm],"relaxed_threshold"->If[aa==2,ToString[1-1/n^2,InputForm],Null],"relaxed_q"->If[aa==2,ToString[N[relax,wp-10],InputForm],Null],"relaxed_increment"->If[aa==2,ToString[N[relax-q,wp-10],InputForm],Null],"relaxed_increment_scaled"->If[aa==2,ToString[N[(relax-q)n(r+1),45],InputForm],Null]|>],{wp,{60,120}}];
 <|"N"->n,"A"->aa,"rows"->rows|>]];
f82[t_]:=t(1+Log[t])^(nn-1);
sym82=<|
"inverse_Phi_quadratic"->TrueQ[FullSimplify[Coefficient[Normal[Series[Log[f82[1+v/nn+(2nn-1)v^2/(2nn^3)]],{v,0,2}]],v,2],Assumptions->nn>=2]==0],
"inverse_Phi_linear"->TrueQ[FullSimplify[Coefficient[Normal[Series[Log[f82[1+v/nn+(2nn-1)v^2/(2nn^3)]],{v,0,2}]],v,1],Assumptions->nn>=2]==1],
"Gamma_density_first_correction"->TrueQ[FullSimplify[Coefficient[Normal[Series[(1/h^2-1)Log[1+x h]-x/h+x^2/2-h^2/12,{h,0,1}]],h,1]-(x^3/3-x)]===0],
"Gamma_CDF_correction_derivative"->TrueQ[FullSimplify[D[(1-x^2)Exp[-x^2/2]/3,x]-(x^3/3-x)Exp[-x^2/2]]==0],
"Fenchel_midpoint_intercept"->TrueQ[FullSimplify[t-f82[t]/D[f82[t],t]==t(nn-1)/(nn+Log[t]),Assumptions->nn>=2&&t>1]],
"inverse_threshold_derivative"->TrueQ[FullSimplify[f82[t]/D[f82[t],t]==t(1+Log[t])/(nn+Log[t]),Assumptions->nn>=2&&t>1]],
"Gamma_low_piece"->TrueQ[FullSimplify[Integrate[Exp[s-q] Exp[-s]s^(nn-1)/Gamma[nn],{s,0,q},Assumptions->nn>=2&&q>0]==Exp[-q]q^nn/Gamma[nn+1],Assumptions->nn>=2&&q>0]]
|>;
dec82[s_String]:=Module[{p=StringSplit[StringReplace[s,"E"->"e"],"e"],v,ss,ee},ss=StringSplit[p[[1]],"."];v=If[Length[ss]==1,ToExpression[ss[[1]],InputForm],ToExpression[StringJoin[ss],InputForm]/10^StringLength[ss[[2]]]];ee=If[Length[p]>1,ToExpression[p[[2]],InputForm],0];v 10^ee];
cert82[z_]:=Module[{n=z["N"],aa=z["A"],q=dec82[z["candidate"]],lo,hi,v}, {lo,hi}=dec82/@z["mean_interval"];v=1+mean82[n,aa,1,N[q,90],90];<|"N"->n,"A"->aa,"side"->z["side"],"candidate"->z["candidate"],"mean"->ToString[N[v,60],InputForm],"inside_directed_interval"->TrueQ[N[lo,90]<v<N[hi,90]],"correct_sign"->TrueQ[If[z["side"]=="upper",v<1,v>1]]|>];
src82=Import[input82,"RawJSON"];certsrc82=Import[FileNameJoin[{DirectoryName[input82],"certificates82.json"}],"RawJSON"];
LaunchKernels[12];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[mean82,run82,dec82,cert82];tt82=AbsoluteTime[];
r82=ParallelMap[TimeConstrained[run82[#],240,"TIMEOUT"]&,src82["rows"],Method->"FinestGrained"];c82=ParallelMap[cert82,certsrc82["rows"],Method->"FinestGrained"];CloseKernels[];
pass82=And@@Values[sym82]&&AllTrue[r82,AssociationQ]&&AllTrue[c82,TrueQ[#["inside_directed_interval"]]&&TrueQ[#["correct_sign"]]&];
Export[FileNameJoin[{out82,"verification82.json"}],<|"version"->$Version,"numerical_strings"->"Working-precision outputs, not certified decimal digits; acceptance uses 60/120 stability and separate directed root brackets","workers"->12,"kernels_closed"->True,"seconds"->N[AbsoluteTime[]-tt82],"pass"->TrueQ[pass82],"symbolic"->sym82,"rows"->r82,"certificate_checks"->c82|>,"RawJSON"];
Print[If[TrueQ[pass82],"PASS ","FAIL "],Length[r82]," optimized corner roots, each at60/120digits"];
If[!TrueQ[pass82],Print[sym82];Print[Cases[r82,_String]];Exit[3]];
