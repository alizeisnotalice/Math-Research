(* r081: independent Wolfram verification; all input rationals remain exact. *)
a81=Rest[$ScriptCommandLine];in81=a81[[1]];out81=a81[[2]];If[!DirectoryQ[out81],CreateDirectory[out81,CreateIntermediateDirectories->True]];
ClearAll[rat81,floor81,ceil81,exp81,cutiv81,chi81,nodes81,cover81,dil81,ext81];
rat81[x_String]:=If[StringMatchQ[x,RegularExpression["-?[0-9]+(/[0-9]+)?"]],ToExpression[x,InputForm],Throw["invalid rational"]];
scale81=10^100;floor81[x_]:=Floor[scale81 x]/scale81;ceil81[x_]:=-floor81[-x];
exp81[x_]:=Module[{b=0,t=x,term=1,lo=1,hi,nxt},While[t>1/2,t=t/2;b++];Do[term=term*t/k;lo+=term,{k,1,80}];nxt=term*t/81;hi=lo+nxt/(1-t/82);{lo,hi}={floor81[lo],ceil81[hi]};Do[{lo,hi}={floor81[lo^2],ceil81[hi^2]},{b}];{lo,hi}];
cutiv81[z_]:=Which[z<=0,{0,0},z>=1,{1,1},z>1/2,Reverse[1-cutiv81[1-z]],True,With[{e=exp81[1/z-1/(1-z)]},{floor81[1/(1+e[[2]])],ceil81[1/(1+e[[1]])]}]];
chi81[z_]:=Which[z<=0,0,z>=1,1,True,1/(1+Exp[1/z-1/(1-z)])];
nodes81[j_]:=Ceiling[Sqrt[2^(j+2)/3]]+1;
cover81[z_]:=Catch[Module[{k=z["k"],jstop,lo=0,hi=0,iv,tail,jf,bs,rows},jstop=k+4;jf=sourcejf81;If[jf!=(rat81/@z["source_modular_interval"]),Throw["actual source modular"]];
Do[iv=cutiv81[2^(k-4-j)];lo+=nodes81[j] iv[[1]];hi+=nodes81[j] iv[[2]],{j,0,jstop-1}];
tail=16nodes81[jstop]/2^256;bs={floor81[jf[[1]]lo],ceil81[jf[[2]](hi+tail)]};If[bs!=(rat81/@z["budget_interval"])||tail!=rat81[z["tail_weighted_node_sum_upper"]],Throw["budget or tail"]];
rows=Table[With[{v=N[(1+5Log[2])/8 Sum[nodes81[j]chi81[2^(k-4-j)],{j,0,jstop-1}],wp]},<|"precision"->wp,"prefix_budget"->ToString[N[v,55],InputForm],"normalized_by_eta_sqrt"->ToString[N[v 2^(-k/2),45],InputForm]|>],{wp,{60,120}}];
<|"k"->k,"exact_pass"->True,"prefix_families"->jstop,"rows"->rows|>]];
dil81[z_]:=Module[{n=z["N"],eta,tau,ze,p,q,lp,lm,tc,rows},eta=1/n^2;tau=-Log[1-eta];ze=Log[1+Log[2]]/Log[2];p=1+ze(n-1);q=1/(p-1);lp=Log[1+eta]+(n-1)Log[1+Log[1+eta]];lm=tau+(n-1)Log[1+tau];tc=n(q-Log[1+q]);rows=Table[Association["precision"->wp,Sequence@@Table[f[[1]]->ToString[N[f[[2]],wp],InputForm],{f,{{"log_dilation_plus",lp},{"log_dilation_minus",lm},{"cap_branch_increment",p tau},{"tau_critical",tc}}}],"cap_branch"->TrueQ[N[tau,wp]<=N[tc,wp]]],{wp,{60,120}}];<|"N"->n,"rows"->rows|>];
ext81[z_]:=Catch[Module[{n=z["N"],tl,th,blo=1,bhi=1,slo=1,shi=1,rb,first,iv,t,eta,coeff,rr,ymax,integrand,def,lv,rows={}},
tl=Sum[2/( (2k+1)3^(2k+1)),{k,0,109}];th=tl+2/(221 3^221)/(1-1/9);{tl,th}={floor81[tl],ceil81[th]};
Do[blo=floor81[blo(n+k)/(k+1)tl^2/((2n+2k+1)(2n+2k+2))];bhi=ceil81[bhi(n+k)/(k+1)th^2/((2n+2k+1)(2n+2k+2))];slo=floor81[slo+blo];shi=ceil81[shi+bhi],{k,0,39}];rb=th^2/(2(2n+1));first=ceil81[bhi(n+40)/41 th^2/((2n+81)(2n+82))];iv={slo,ceil81[shi+first/(1-rb)]};If[iv!=(rat81/@z["normalized_volume_interval"]),Throw["positive volume certificate"]];
Do[t=N[Log[2],wp];eta=N[1/n^2,wp];coeff=N[FoldList[#1 (n+#2)/(#2+1) Log[2]^2/((2n+2#2+1)(2n+2#2+2))&,1,Range[0,58]],wp];rr=Total[coeff];ymax=-2n Log[1-Log[1+eta]/t];
integrand[v_?NumericQ]:=Module[{y=ymax v,zz,density,cut},zz=(Exp[t(1-Exp[-y/(2n)])]-1)/eta;cut=chi81[zz];density=Sum[coeff[[k+1]](1+k/n)Exp[-(1+k/n)y],{k,0,59}]/rr;ymax(1-cut)density];
def=NIntegrate[integrand[v],{v,0,1/4,1/2,3/4,1},WorkingPrecision->wp,AccuracyGoal->Floor[wp/2],PrecisionGoal->Floor[wp/2],Method->{"GlobalAdaptive","SymbolicProcessing"->0},MaxRecursion->20];
lv=2n Log[t]-LogGamma[2n+1]+Log[rr];AppendTo[rows,<|"precision"->wp,"normalized_volume_series"->ToString[N[rr,55],InputForm],"log_exterior_volume"->ToString[N[lv,55],InputForm],"log_exterior_ratio_upper"->ToString[N[lv+(n-1)t,55],InputForm],"smooth_to_hard_volume"->ToString[N[1-def,55],InputForm],"scaled_smoothing_deficit"->ToString[N[n def,55],InputForm],"log_modular_relative_error_upper"->ToString[N[n Log[2t]-LogGamma[n],55],InputForm]|>],{wp,{60,120}}];
<|"N"->n,"exact_pass"->True,"rows"->rows|>]];
l2lo81=Sum[2/((2k+1)3^(2k+1)),{k,0,109}];l2hi81=l2lo81+2/(221 3^221)/(1-1/9);sourcejf81={floor81[(1+5l2lo81)/8],ceil81[(1+5l2hi81)/8]};
soft81=Import[FileNameJoin[{in81,"soft81.json"}],"RawJSON"];exterior81=Import[FileNameJoin[{in81,"exterior81.json"}],"RawJSON"];
expr81=1/(1+Exp[1/z-1/(1-z)]);
sym81=<|
"shared_exact_cutoff_intervals"->And@@Table[cutiv81[2^(-z["m"]) ]==(rat81/@z["interval"]),{z,soft81["flat_weight_intervals"]}],
"cutoff_derivative"->TrueQ[FullSimplify[D[expr81,z]==(1/z^2+1/(1-z)^2)expr81(1-expr81),Assumptions->0<z<1]],
"cutoff_symmetry"->TrueQ[FullSimplify[expr81+(expr81/.z->1-z)==1,Assumptions->0<z<1]],
"endpoint_limits_0_to_3"->And@@Table[TrueQ[TimeConstrained[Limit[D[expr81,{z,k}],z->0,Direction->"FromAbove"],15,$Failed]===0],{k,0,3}],
"dyadic_node_growth"->TrueQ[Resolve[ForAll[x,x>0,Ceiling[Sqrt[2]x]+1<=2(Ceiling[x]+1)],Reals]],
"geometric_volume_ratio"->TrueQ[Resolve[ForAll[{n,k},n>=2&&k>=0,(n+k)/((k+1)(2n+2k+1)(2n+2k+2))<=1/(2(2n+1))],Reals]],
"exterior_volume_two"->TrueQ[FullSimplify[LaplaceTransform[1-Cosh[t]+t Sinh[t]/2,t,s,Assumptions->s>1]==1/(s(s^2-1)^2)]],
"exterior_volume_two_log2"->TrueQ[FullSimplify[1-Cosh[Log[2]]+Log[2]Sinh[Log[2]]/2==(3Log[2]-2)/8]],
"r080_cover_margin"->TrueQ[229/768==7/12-73/256],
"tail_e_bounds"->TrueQ[(rat81/@soft81["e_interval"])=={Sum[1/i!,{i,0,10}],Sum[1/i!,{i,0,10}]+1/11!/(1-1/12)}&&2<rat81[soft81["e_interval"][[1]]]<rat81[soft81["e_interval"][[2]]]<11/4&&(11/4)^2<8]
|>;
(* Ceiling inequality is elementary from ceil(sqrt(2)x)<=2ceil(x); avoid unsupported quantifier elimination. *)
If[!TrueQ[sym81["dyadic_node_growth"]],sym81["dyadic_node_growth"]=TrueQ[Resolve[ForAll[{x,m},m>=1&&0<x<=m,Sqrt[2]x<=2m],Reals]]];
LaunchKernels[12];ParallelEvaluate[$HistoryLength=0];DistributeDefinitions[rat81,sourcejf81,scale81,floor81,ceil81,exp81,cutiv81,chi81,nodes81,cover81,dil81,ext81];tt81=AbsoluteTime[];
c81=ParallelMap[TimeConstrained[cover81[#],180,"TIMEOUT"]&,soft81["cover"],Method->"FinestGrained"];
d81=ParallelMap[dil81,soft81["dilation"],Method->"FinestGrained"];
e81=ParallelMap[TimeConstrained[ext81[#],240,"TIMEOUT"]&,exterior81["rows"],Method->"FinestGrained"];
CloseKernels[];
pass81=And@@Values[sym81]&&AllTrue[Join[c81,e81],AssociationQ[#]&&TrueQ[#["exact_pass"]]&]&&AllTrue[d81,AssociationQ];
Export[FileNameJoin[{out81,"verification81.json"}],<|"version"->$Version,"workers"->12,"kernels_closed"->True,"seconds"->N[AbsoluteTime[]-tt81],"symbolic"->sym81,"pass"->TrueQ[pass81],"cover"->c81,"dilation"->d81,"exterior"->e81|>,"RawJSON"];
Print[If[TrueQ[pass81],"PASS ","FAIL "],Length[c81]," flat cover scales; ",Length[d81]," dimensions; ",Length[e81]," exterior integrals"];
If[!TrueQ[pass81],Print[sym81];Print[Cases[Join[c81,e81],_String]];Exit[3]];
