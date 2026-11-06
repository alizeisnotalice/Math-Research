(* r041: unequal Gamma--binomial blocks n1+n2=64.  Source and receiver
   kernels are rebuilt for each block; all 65 anchors are tested separately.
   Floating split-quadrature evidence only.
   wolframscript -file code/unbalanced_blocks41.wl [out.json] [seed.json] *)
root41=DirectoryName[$InputFileName];args41=$ScriptCommandLine;repo41=DirectoryName[root41];seed41=If[Length[args41]>=3,args41[[3]],FileNameJoin[{repo41,"rounds","r030","receiver_seed.json"}]];out41=If[Length[args41]>=2,args41[[2]],FileNameJoin[{$TemporaryDirectory,"unbalanced_blocks41.json"}]];
$r032LibraryMode=True;Block[{$ScriptCommandLine={"unbalanced_blocks41.wl",seed41}},Get[FileNameJoin[{repo41,"code","receiver_continuous32.wl"}]]];$r032LibraryMode=False;

(* r039 source weights and per-anchor ratio template.  No `polyC` or global
   b=32 is reused below: block sizes enter every density and CDF explicitly. *)
cBase41=c; h41=.15; eta41=.2; scale41=1.020; qv41=vv/16.;
wRaw41=w0; w41=wRaw41 Exp[eta41 qv41^2] Total[wRaw41]/Total[wRaw41 Exp[eta41 qv41^2]];
csG41=cBase41+h41 qv41^2; cAt41=Join[{csG41[[1]]},Rest[csG41],Rest[csG41]];
bin41[bb_Integer]:=N[Binomial[bb,#]&/@Range[0,bb]];
polyB41[t_?NumericQ,bb_Integer,cj_?NumericQ,bi_List]:=If[t<=0.,0.,Total[bi Clip[t-cj Range[0,bb],{0.,Infinity}]^(bb-1)]/(Gamma[bb](1.+Exp[-cj])^bb)];
cdfB41[z_?NumericQ,bb_Integer,cj_?NumericQ,bi_List]:=Module[{ww=bi Exp[-cj Range[0,bb]]/(1.+Exp[-cj])^bb},Total[ww GammaRegularized[bb,0.,Clip[z-cj Range[0,bb],{0.,Infinity}]]]];
cone41[anchors_]:=Module[{s=SortBy[anchors,First],xx,best,tot=0.,i},xx=DeleteDuplicates[s[[All,1]]];
 Do[best=Min[Select[s,#[[1]]<=xx[[i]]&][[All,2]]];tot+=(Exp[-xx[[i]]]-If[i<Length[xx],Exp[-xx[[i+1]]],0.])Exp[-best],{i,Length[xx]}];tot];

run41[name_String,n1_Integer,n2_Integer,q_?NumericQ,no_Integer:64,ord_Integer:32]:=Module[
 {dnew,u41,lo41,hi41,bi1,bi2,sourceW,pts,idx,amp,like,bd,pv,cp,sp,tx,tw,aa,roots,z,a,cost=0.,av,ids,an,vol,cm,pmax,mzMax,vzMax},
 dnew=dd+q qv41^2;
 u41=Join[{{0.,0.}},Transpose[{Rest[dnew]+Rest[vv],Rest[dnew]-Rest[vv]}],Transpose[{Rest[dnew]-Rest[vv],Rest[dnew]+Rest[vv]}]];
 lo41=Min[u41];cm=Max[cAt41];pmax=1./(1.+Exp[cm]);mzMax=1.+cm pmax;vzMax=1.+cm^2 pmax(1.-pmax);
 hi41=Max[u41]+Max[n1,n2] mzMax+16. Sqrt[Max[n1,n2] vzMax];
 bi1=bin41[n1];bi2=bin41[n2];sourceW=scale41 Join[{w41[[1]]},Rest[w41],Rest[w41]];
 pts=Range[2ng-1];idx[k_Integer]:=If[k<=ng,k,k-ng+1];
 amp[t1_?NumericQ]:=Table[sourceW[[k]] polyB41[t1-u41[[k,1]],n1,cAt41[[k]],bi1] Exp[u41[[k,1]]+u41[[k,2]]-tau],{k,pts}];
 like[t2_?NumericQ,avec_]:=Total[avec Table[polyB41[t2-u41[[k,2]],n2,cAt41[[k]],bi2],{k,pts}]];
 bd[lev_?NumericQ,avec_]:=Module[{l=lo41,r=hi41,m},If[like[l,avec]>=lev,Return[l]];If[like[r,avec]<=lev,Return[r]];
   Do[m=(l+r)/2.;If[like[m,avec]<lev,l=m,r=m],{45}];(l+r)/2.];
 pv[t1_?NumericQ,t2_?NumericQ]:=Table[polyB41[t1-u41[[k,1]],n1,cAt41[[k]],bi1] polyB41[t2-u41[[k,2]],n2,cAt41[[k]],bi2] Exp[u41[[k,1]]+u41[[k,2]]-t1-t2],{k,pts}];
 cp[t1_?NumericQ,l_?NumericQ,r_?NumericQ,hh_?NumericQ]:=Module[{ccost=0.,aav=ConstantArray[0.,2ng-1],fac,k},If[r<=l,Return[{ccost,aav}]];
   ccost=Exp[tau-t1](Exp[-l]-Exp[-r])phi[hh];Do[fac=polyB41[t1-u41[[k,1]],n1,cAt41[[k]],bi1] Exp[u41[[k,1]]-t1];
    aav[[k]]=fac(cdfB41[r-u41[[k,2]],n2,cAt41[[k]],bi2]-cdfB41[l-u41[[k,2]],n2,cAt41[[k]],bi2]),{k,pts}];{ccost,aav}];
 sp[t1_?NumericQ,l_?NumericQ,r_?NumericQ,avec_,oo_Integer]:=Module[{x,wt,hh,ccost=0.,aav=ConstantArray[0.,2ng-1]},If[r<=l,Return[{ccost,aav}]];{x,wt}=quad[oo,l,r];
   Do[hh=profile[like[x[[i]],avec]];ccost+=wt[[i]] Exp[tau-t1-x[[i]]]phi[hh];aav+=wt[[i]]hh pv[t1,x[[i]]],{i,Length[x]}];{ccost,aav}];
 {tx,tw}=Module[{pan=Ceiling[no/16.],edge,x,wt,outx={},outw={}},edge=Subdivide[lo41,hi41,pan];Do[{x,wt}=quad[16,edge[[i]],edge[[i+1]]];outx=Join[outx,x];outw=Join[outw,wt],{i,pan}];{outx,outw}];av=ConstantArray[0.,2ng-1];
 Do[aa=amp[tx[[i]]];roots=bd[#,aa]&/@{1.,N[n],phiPrime};{z,a}=cp[tx[[i]],roots[[1]],roots[[2]],1.];cost+=tw[[i]]z;av+=tw[[i]]a;
   {z,a}=cp[tx[[i]],roots[[3]],hi41,cap];cost+=tw[[i]]z;av+=tw[[i]]cap a;{z,a}=sp[tx[[i]],roots[[2]],roots[[3]],aa,ord];cost+=tw[[i]]z;av+=tw[[i]]a,{i,Length[tx]}];
 ids=Flatten@Position[av,_?(#>1.00001&)];an=u41[[ids]];vol=If[Length[an]>0,cone41[an],0.];
 <|"name"->name,"blocks"->{n1,n2},"q"->q,"source_scale"->scale41,"eta"->eta41,"h"->h41,"outer_nodes"->no,"soft_order"->ord,
 "scaled_cost"->cost,"actual_logJ"->Log[cost]-tau,"anchor_averages"->av,"min_anchor_average"->Min[av],"accepted_anchor_count"->Length[ids],"accepted_anchor_indices"->ids,
 "accepted_cone_volume"->vol,"accepted_cone_log_score"->If[vol>0,Log[vol]-(Log[cost]-tau),-Infinity]|>];

If[!TrueQ[$unbalanced41LibraryMode],
cases41={{"balanced_q0",32,32,0.},{"balanced_qminus",32,32,-.08},{"b24_40_q0",24,40,0.},{"b24_40_qminus",24,40,-.08},{"b16_48_q0",16,48,0.},{"b16_48_qminus",16,48,-.08}};
initial41=(run41[#[[1]],#[[2]],#[[3]],#[[4]]]& /@ cases41);
top41=Take[Ordering[initial41[[All,"accepted_cone_log_score"]],-2],2];
valid41=Table[Join[<|"selected_by_N64_score"->initial41[[k,"name"]]|>,run41[initial41[[k,"name"]],initial41[[k,"blocks",1]],initial41[[k,"blocks",2]],initial41[[k,"q"]],128,32]],{k,top41}];
res41=<|"status"->"r041 unequal-block source/receiver reconstruction; floating direct split quadrature, not a certificate",
 "wolfram"->$Version,"n"->n,"source_weight_formula"->"w_j exp(.2(v_j/16)^2) sum(w)/sum(w exp(.2(v_j/16)^2))","receiver_template"->"c_j=c0+.15(v_j/16)^2",
 "geometry_formula"->"d_j -> d_j + q(v_j/16)^2 for q=0 or -.08, applied to source shifts and receiver anchors",
 "anchor_accounting"->"all 65 oriented anchors evaluated separately; no mirror averaging or division by 2", "initial_N64_soft32"->initial41,"selected_N128_soft32"->valid41,
 "scalar_reference_log_score"->71.71499208794934,"selection"->"actual Phi64 plus finite accepted-cone score; individual anchor average >1+1e-5"|>;
Export[out41,res41,"RawJSON"];Print[<|"output"->out41,"cases"->Length[initial41],"refinements"->Length[valid41]|>];
];
