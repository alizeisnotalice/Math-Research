(* r050: gauge-fixed relative receiver geometry probe, n=64.
   Every trial keeps u_1=(0,0); source H is rebuilt after the relative
   geometry change.  Floating split quadrature only. *)
root50=DirectoryName[$InputFileName];repo50=DirectoryName[root50];out50=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/relative50.json"];
$r032LibraryMode=True;Block[{$ScriptCommandLine={"receiver_continuous32.wl"}},Get[FileNameJoin[{repo50,"code","receiver_continuous32.wl"}]]];$r032LibraryMode=False;
cBase50=c; q50=vv/16.; srcC50=cBase50+.15q50^2; recUniform50=ConstantArray[cBase50+.10,ng];
wRaw50=w0;w50=wRaw50 Exp[.2q50^2] Total[wRaw50]/Total[wRaw50 Exp[.2q50^2]];srcScale50=1.020;
setGeom50[span_?NumericQ,curv_?NumericQ]:=Module[{dnew=dd+curv q50^2,vnew=(1+span)vv,cm,pp},
 u=Join[{{0.,0.}},Transpose[{Rest[dnew]+Rest[vnew],Rest[dnew]-Rest[vnew]}],Transpose[{Rest[dnew]-Rest[vnew],Rest[dnew]+Rest[vnew]}]];
 cm=Max[srcC50];pp=1/(1+Exp[cm]);lo=Min[u];hi=Max[u]+b(1+cm pp)+16 Sqrt[b(1+cm^2 pp(1-pp))];];
pts50[j_Integer]:=If[j==1,{u[[1]]},{u[[j]],u[[j+ng-1]]}];
mixAmp50[t1_?NumericQ]:=Module[{aa=ConstantArray[0.,2ng-1],z},
 aa[[1]]=srcScale50 w50[[1]] polyC[t1-u[[1,1]],srcC50[[1]]] Exp[Total[u[[1]]]-tau];
 Do[z=u[[j]];aa[[j]]=srcScale50 w50[[j]] polyC[t1-z[[1]],srcC50[[j]]] Exp[Total[z]-tau];z=u[[j+ng-1]];aa[[j+ng-1]]=srcScale50 w50[[j]] polyC[t1-z[[1]],srcC50[[j]]] Exp[Total[z]-tau],{j,2,ng}];aa];
like50[t2_?NumericQ,aa_]:=Total[aa Table[polyC[t2-u[[k,2]],srcC50[[If[k<=ng,k,k-ng+1]]]],{k,2ng-1}]];
bd50[lev_?NumericQ,aa_]:=Module[{l=lo,r=hi,m},If[like50[l,aa]>=lev,Return[l]];If[like50[r,aa]<=lev,Return[r]];Do[m=(l+r)/2.;If[like50[m,aa]<lev,l=m,r=m],{45}];(l+r)/2.];
pv50[t1_?NumericQ,t2_?NumericQ,rc_]:=Module[{av=ConstantArray[0.,ng],z},
 Do[Do[av[[j]]+=polyC[t1-z[[1]],rc[[j]]] polyC[t2-z[[2]],rc[[j]]] Exp[Total[z]-t1-t2],{z,pts50[j]}],{j,ng}];av];
cp50[t1_?NumericQ,l_?NumericQ,r_?NumericQ,hh_?NumericQ,rc_]:=Module[{cost=0.,av=ConstantArray[0.,ng],z,fac},If[r<=l,Return[{cost,av}]];cost=Exp[tau-t1](Exp[-l]-Exp[-r])phi[hh];
 Do[Do[fac=polyC[t1-z[[1]],rc[[j]]] Exp[z[[1]]-t1];av[[j]]+=fac(gbcdfC[r-z[[2]],rc[[j]]]-gbcdfC[l-z[[2]],rc[[j]]]),{z,pts50[j]}],{j,ng}];{cost,av}];
sp50[t1_?NumericQ,l_?NumericQ,r_?NumericQ,aa_,rc_,ord_Integer]:=Module[{x,wt,hh,cost=0.,av=ConstantArray[0.,ng]},If[r<=l,Return[{cost,av}]];{x,wt}=quad[ord,l,r];Do[hh=profile[like50[x[[i]],aa]];cost+=wt[[i]]Exp[tau-t1-x[[i]]]phi[hh];av+=wt[[i]]hh pv50[t1,x[[i]],rc],{i,Length[x]}];{cost,av}];
cone50[anchors_]:=Module[{s=SortBy[anchors,First],xx,best,tot=0.,i},xx=DeleteDuplicates[s[[All,1]]];Do[best=Min[Select[s,#[[1]]<=xx[[i]]&][[All,2]]];tot+=(Exp[-xx[[i]]]-If[i<Length[xx],Exp[-xx[[i+1]]],0.])Exp[-best],{i,Length[xx]}];tot];
run50[name_String,span_?NumericQ,curv_?NumericQ,recName_String,rc_,no_Integer:128,ord_Integer:32]:=Module[{tx,tw,aa,roots,z,a,cost=0.,av=ConstantArray[0.,ng],ids,an,vol},setGeom50[span,curv];{tx,tw}=outerQuad[no];
 Do[aa=mixAmp50[tx[[i]]];roots=bd50[#,aa]&/@{1.,N[n],phiPrime};{z,a}=cp50[tx[[i]],roots[[1]],roots[[2]],1.,rc];cost+=tw[[i]]z;av+=tw[[i]]a;{z,a}=cp50[tx[[i]],roots[[3]],hi,cap,rc];cost+=tw[[i]]z;av+=tw[[i]]cap a;{z,a}=sp50[tx[[i]],roots[[2]],roots[[3]],aa,rc,ord];cost+=tw[[i]]z;av+=tw[[i]]a,{i,Length[tx]}];
 ids=Flatten@Position[av/coeff,_?(#>1.00001&)];an=Flatten[Table[pts50[j],{j,ids}],1];vol=If[Length[an]>0,cone50[an],0.];<|"geometry"->name,"span"->span,"curvature"->curv,"receiver_parameter"->recName,"outer_nodes"->no,"soft_order"->ord,"finite_domain"->{lo,hi},"source_cap"->cap,"scaled_cost"->cost,"actual_logJ"->Log[cost]-tau,"averages_by_group"->av/coeff,"accepted_group_ids"->ids,"min_group_average"->Min[av/coeff],"strict_margin"->Min[av/coeff]-1,"accepted_groups"->Length[ids],"accepted_nodes"->Length[an],"all_65_nodes_pass"->(Length[ids]==ng),"accepted_cone_volume"->vol,"accepted_cone_log_score"->If[vol>0,Log[vol]-(Log[cost]-tau),-Infinity]|>];
tailCache50=<||>;
tailDiagnostics50[row_]:=Module[{key=ToString[Lookup[row,{"span","curvature","receiver_parameter"}],InputForm],rc,bcost,bavg,tails,upper,cc,z0,z,wh,wl,ta},
 If[KeyExistsQ[tailCache50,key],Return[Join[row,tailCache50[key]]]];setGeom50[row["span"],row["curvature"]];rc=If[row["receiver_parameter"]=="matched_cj",srcC50,recUniform50];
 wh=Rationalize[hi,0];wl=Rationalize[lo,0];ta=Rationalize[tau,0];
 upper[zz_,cj_]:=Module[{cc=Rationalize[cj,0],ww},ww=N[(Binomial[b,#]Exp[-cc #]/(1+Exp[-cc])^b)&/@Range[0,b],60];Total[ww GammaRegularized[b,Clip[N[zz-cc Range[0,b],60],{0,Infinity}]]]];
 tails=Table[z0=First[pts50[j]];z=Rationalize[z0,0];upper[wh-z[[1]],rc[[j]]]+upper[wh-z[[2]],rc[[j]]],{j,ng}];
 bavg=cap Max[tails];bcost=N[Exp[ta]cap(1+Log[cap])^(n-1)(2Exp[-wl-wh]-Exp[-2wh]),60];
 tailCache50[key]=<|"scaled_budget_tail_upper_formula_value"->N[bcost,25],"omitted_average_upper_formula_value"->N[bavg,25],"tail_evaluation_precision"->60,"tail_scope"->"analytic cap/union probability bounds evaluated numerically with upper Gamma tails; not a total quadrature/rounding certificate"|>;Join[row,tailCache50[key]]];
If[!TrueQ[$relative50LibraryMode],
cases50={{"baseline",0.,0.},{"angular_expand",.10,0.},{"angular_contract",-.10,0.},{"outer_curve_in",0.,-.08}};
rows50=Flatten[Table[{run50[ca[[1]],ca[[2]],ca[[3]],"matched_cj",srcC50],run50[ca[[1]],ca[[2]],ca[[3]],"uniform_c0_plus_.10",recUniform50]},{ca,cases50}],1];
best50=First@SortBy[Select[rows50,#"receiver_parameter"=="matched_cj"&],-#"accepted_cone_log_score"&];ref50=run50[best50["geometry"],best50["span"],best50["curvature"],"matched_cj",srcC50,256,64];
bestAny50=First@SortBy[rows50,-#"accepted_cone_log_score"&];refAny50=run50[bestAny50["geometry"],bestAny50["span"],bestAny50["curvature"],bestAny50["receiver_parameter"],If[bestAny50["receiver_parameter"]=="matched_cj",srcC50,recUniform50],256,64];
rows50=tailDiagnostics50/@rows50;ref50=tailDiagnostics50[ref50];refAny50=tailDiagnostics50[refAny50];
res50=<|"wolfram"->$Version,"status"->"r050 gauge-fixed relative receiver geometry probe; floating continuous branch-split quadrature, not a certificate","n"->n,"gauge_fix"->"u_1=(0,0); no common translation parameter; all u changes are relative span/curvature deformations","source"->"r039 weights exp(.2(v/16)^2)-normalized, scale 1.020, c_source_j=c0+.15(v/16)^2 rebuilt at each geometry","trials_N128_soft32"->rows50,"best_matched_N256_soft64"->ref50,"best_overall_N256_soft64"->refAny50,"scalar_reference_log_score"->71.71499208794934,"scope"->"For each geometry matched_cj and uniform receiver-c control share exactly the same rebuilt source H; cost differences between those two are numerical integration noise only. Full halo is represented only by accepted finite Pareto cones."|>;
Export[out50,res50,"RawJSON"];Print[<|"rows"->rows50[[All,{"geometry","receiver_parameter","accepted_cone_log_score","strict_margin","all_65_nodes_pass"}]],"ref"->ref50|>];

];
