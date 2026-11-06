(* Numerical test of the clipped exact quadratic-Poisson potential.
   All statements from this script are floating diagnostics, not certificates. *)
repo="/Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal/strong-maximal-endpoint";
src=Import[FileNameJoin[{repo,"code","payoff_bellman33.wl"}],"Text"];
ToExpression[StringTake[src,First[First[StringPosition[src,"valueAt[values_List"]]]-1]];
mm=1+ProductLog[Exp[-1]];cc=N[mm,50];pp=1/(1+Exp[cc]);
Clear[stats,make];
stats[nn_Integer,tt_Integer,rin_?NumericQ]:=stats[nn,tt,rin]=Module[{rr=SetPrecision[rin,50],ks=Range[0,tt],x,bb,ww,s0,s1,s2,den,ff,fp,fc,fcp,rs,qq,call,vv,dd,uu,up,rho,xx,bc,bp},
 If[rr<=0,Return[N[{1-rr/nn+tt cc/nn,0,-1/nn,1-rr/nn+tt cc/nn,-1/nn,0,0,0},17]]];
 x=Clip[rr-cc ks,{0,Infinity}];bb=N[Binomial[tt,#]&/@ks,50];ww=bb Exp[-cc ks]/(1+Exp[-cc])^tt;
 den=(1+Exp[-cc])^tt Gamma[tt];s0=Total[bb x^(tt-1)];s1=Total[bb x^(tt-2)];
 ff=Exp[-rr]s0/den;fp=Exp[-rr]((tt-1)s1-s0)/den;
 fc=Exp[-rr]Total[bb(x^(tt+1)/(tt+1)+cc^2 ks(tt-ks)x^(tt-1)/(tt-1))]/den;
 fcp=Exp[-rr]Total[bb(x^tt-x^(tt+1)/(tt+1)+cc^2 ks(tt-ks)((tt-1)x^(tt-2)-x^(tt-1))/(tt-1))]/den;
 rs=Exp[-rr]Total[bb x^tt]/(den tt);qq=Total[ww GammaRegularized[tt,x]];
 call=Total[ww(tt GammaRegularized[tt+1,x]-(rr-cc ks)GammaRegularized[tt,x])];
 vv=rs+qq+call/nn;dd=rs+qq/nn;
 uu=pp/8 (ff(rr-tt cc)^2/tt+(HarmonicNumber[tt]-1)fc/tt);
 up=pp/8 ((fp(rr-tt cc)^2+2ff(rr-tt cc))/tt+(HarmonicNumber[tt]-1)fcp/tt);
 rho=nn+ProductLog[-nn Exp[-nn]];xx=rr-rho;
 bc=If[xx<=0,0,Exp[-rr]xx^tt/Gamma[tt+1]+(rr/nn-1)GammaRegularized[tt,0,xx]-tt/nn GammaRegularized[tt+1,0,xx]];
 bp=If[xx<=0,0,GammaRegularized[tt,0,xx]/nn-Exp[-rr]xx^tt/Gamma[tt+1]];
 N[{vv,uu,-dd,1-rr/nn+tt cc/nn+bc,-1/nn+bp,ff,up,fp},17]];
make[nn_,tt_,alpha_,beta_,hh_]:=Module[{grid,st,prev,cur,iv,pv,ans,active,viol,res,idx,core},
 grid=N[Range[0,Ceiling[(tt cc+5Sqrt[tt cc])/hh]]hh];
 st=stats[nn,tt,#]&/@grid;prev=stats[nn,tt-1,#]&/@grid;
 iv=st[[All,1]]+alpha st[[All,2]]+beta N[HarmonicNumber[tt,2]]st[[All,6]];
 pv=prev[[All,1]]+alpha prev[[All,2]]+beta N[HarmonicNumber[tt-1,2]]prev[[All,6]];
 cur=MapThread[Min,{iv,st[[All,4]]}];pv=MapThread[Min,{pv,prev[[All,4]]}];
 active=Flatten@Position[MapThread[Less,{iv,st[[All,4]]}],True];
 viol=Max[0,Max[Differences[pv]],Max[Differences[cur]]];
 ans=compiledStage[pv,nn,tt,hh,N[mm]];res=Take[ans,Length[cur]]-cur;
 idx=First@Ordering[res,-1];core=Flatten@Position[(#<=2Sqrt[tt Log[tt]]&)/@Abs[grid-tt N[mm]],True];
 <|"N"->nn,"t"->tt,"alpha"->alpha,"beta"->beta,"h"->hh,"grid_monotonicity_positive_increment"->viol,
 "max_residual"->res[[idx]],"at_r"->grid[[idx]],"residual_over_f"->res[[idx]]/st[[idx,6]],
 "core_max_residual"->Max[res[[core]]],"active_points"->Length[active],
 "max_active_inner_derivative"->Max[st[[active,3]]+alpha st[[active,7]]+beta N[HarmonicNumber[tt,2]]st[[active,8]]],
 "note"->"T uses PWL future interpolation and r033 global root; only reliable as diagnostic when monotonicity check passes."|>];
rows=Flatten[Table[make[256,tt,1.1,bb,hh],{tt,{16,64}},{bb,{.1,1.}},{hh,{.05,.025}}]];
Export["/Users/zhengzhihao/Desktop/T/work/r035/patch_probe35_error.json",<|"wolfram"->$Version,"rows"->rows|>,"RawJSON"];
Print[rows];
