(* r040: synchronous source/receiver curve geometry deformation, n=64.
   Floating split-quadrature evidence only.  Run:
   wolframscript -file lower_geometry40.wl [out.json] [seed.json] *)
root40=DirectoryName[$InputFileName];args40=$ScriptCommandLine;repo40=DirectoryName[root40];seed40=If[Length[args40]>=3,args40[[3]],FileNameJoin[{repo40,"rounds","r030","receiver_seed.json"}]];out40=If[Length[args40]>=2,args40[[2]],FileNameJoin[{$TemporaryDirectory,"lower_geometry40.json"}]];
$r032LibraryMode=True;Block[{$ScriptCommandLine={"lower_geometry40.wl",seed40}},Get[FileNameJoin[{repo40,"code","receiver_continuous32.wl"}]]];$r032LibraryMode=False;

(* The fixed r039 input parameters.  The eta normalization is exact as a
   formula (the following evaluation is floating, like the quadrature). *)
cBase40=c; h40=.15; eta40=.2; scale40=1.020;
wRaw40=w0; totalW40=Total[wRaw40]; qv40=vv/16.;
w0=wRaw40 Exp[eta40 qv40^2] totalW40/Total[wRaw40 Exp[eta40 qv40^2]];
cs40=cBase40+h40 qv40^2;

(* `setGeometry` changes the actual source shifts as well as the receiver
   anchors.  Coordinates are mirrored after the new d_j is formed. *)
setGeometry[power_Integer,amp_?NumericQ]:=Module[{dnew,cm,pmax,vmax},
 dnew=dd+amp qv40^power;
 u=Join[{{0.,0.}},Transpose[{Rest[dnew]+Rest[vv],Rest[dnew]-Rest[vv]}],
   Transpose[{Rest[dnew]-Rest[vv],Rest[dnew]+Rest[vv]}]];
 cm=Max[cs40]; pmax=1./(1.+Exp[cm]);
 mz=1.+cm pmax; vz=1.+cm^2 pmax(1.-pmax);
 lo=Min[u]; hi=Max[u]+b mz+16. Sqrt[b vz];
];
pts40[j_Integer]:=If[j==1,{u[[1]]},{u[[j]],u[[j+ng-1]]}];
mixAmp40[t1_?NumericQ]:=Module[{aa=ConstantArray[0.,2ng-1],z},
 Do[Do[aa[[If[j==1,1,If[z===pts40[j][[1]],j,j+ng-1]]]]=
     scale40 w0[[j]] polyC[t1-z[[1]],cs40[[j]]] Exp[z[[1]]+z[[2]]-tau],
   {z,pts40[j]}],{j,ng}]; aa];
like40[t2_?NumericQ,aa_]:=Module[{s=aa[[1]] polyC[t2-u[[1,2]],cs40[[1]]],j},
 Do[s+=aa[[j]] polyC[t2-u[[j,2]],cs40[[j]]]+aa[[j+ng-1]] polyC[t2-u[[j+ng-1,2]],cs40[[j]]],{j,2,ng}];s];
bd40[lev_?NumericQ,aa_]:=Module[{l=lo,r=hi,m},If[like40[l,aa]>=lev,Return[l]];If[like40[r,aa]<=lev,Return[r]];
 Do[m=(l+r)/2.;If[like40[m,aa]<lev,l=m,r=m],{45}];(l+r)/2.];
pv40[t1_?NumericQ,t2_?NumericQ]:=Module[{av=ConstantArray[0.,ng],z,q1,q2},
 Do[Do[q1=polyC[t1-z[[1]],cs40[[j]]];q2=polyC[t2-z[[2]],cs40[[j]]];
   av[[j]]+=q1 q2 Exp[z[[1]]+z[[2]]-t1-t2],{z,pts40[j]}],{j,ng}];av];
cp40[t1_?NumericQ,l_?NumericQ,r_?NumericQ,hh_?NumericQ]:=Module[{cost=0.,av=ConstantArray[0.,ng],z,fac},
 If[r<=l,Return[{cost,av}]];cost=Exp[tau-t1](Exp[-l]-Exp[-r])phi[hh];
 Do[Do[fac=polyC[t1-z[[1]],cs40[[j]]] Exp[z[[1]]-t1];
   av[[j]]+=fac (gbcdfC[r-z[[2]],cs40[[j]]]-gbcdfC[l-z[[2]],cs40[[j]]]),{z,pts40[j]}],{j,ng}];{cost,av}];
sp40[t1_?NumericQ,l_?NumericQ,r_?NumericQ,aa_,ord_Integer]:=Module[{x,wt,hh,cost=0.,av=ConstantArray[0.,ng]},
 If[r<=l,Return[{cost,av}]];{x,wt}=quad[ord,l,r];
 Do[hh=profile[like40[x[[i]],aa]];cost+=wt[[i]] Exp[tau-t1-x[[i]]] phi[hh];av+=wt[[i]] hh pv40[t1,x[[i]]],{i,Length[x]}];{cost,av}];
cone40[anchors_]:=Module[{s=SortBy[anchors,First],xx,best,tot=0.,i},
 xx=DeleteDuplicates[s[[All,1]]];Do[best=Min[Select[s,#[[1]]<=xx[[i]]&][[All,2]]];
  tot+=(Exp[-xx[[i]]]-If[i<Length[xx],Exp[-xx[[i+1]]],0.]) Exp[-best],{i,Length[xx]}];tot];
run40[name_String,power_Integer,amp_?NumericQ,no_Integer:64,ord_Integer:32]:=Module[{tx,tw,aa,roots,z,a,cost=0.,av=ConstantArray[0.,ng],ids,an,vol},
 setGeometry[power,amp];{tx,tw}=outerQuad[no];
 Do[aa=mixAmp40[tx[[i]]];roots=bd40[#,aa]&/@{1.,N[n],phiPrime};
  {z,a}=cp40[tx[[i]],roots[[1]],roots[[2]],1.];cost+=tw[[i]]z;av+=tw[[i]]a;
  {z,a}=cp40[tx[[i]],roots[[3]],hi,cap];cost+=tw[[i]]z;av+=tw[[i]]cap a;
  {z,a}=sp40[tx[[i]],roots[[2]],roots[[3]],aa,ord];cost+=tw[[i]]z;av+=tw[[i]]a,
 {i,Length[tx]}];ids=Flatten@Position[av/coeff,_?(#>1.00001&)];
 an=Flatten[Table[pts40[j],{j,ids}],1];vol=If[Length[an]>0,cone40[an],0.];
 <|"name"->name,"deformation"->If[amp==0.,"baseline","d_j + q (v_j/16)^power"],"power"->power,"q"->amp,
   "source_scale"->scale40,"eta"->eta40,"receiver_c_template"->"c0+.15(v/16)^2","outer_nodes"->no,"soft_order"->ord,
   "scaled_cost"->cost,"actual_logJ"->Log[cost]-tau,"min_average"->Min[av/coeff],
   "accepted_group_count"->Length[ids],"accepted_groups_margin_1e5"->ids,"accepted_cone_volume"->vol,
   "accepted_cone_log_score"->If[vol>0,Log[vol]-(Log[cost]-tau),-Infinity]|>];

(* q=0 is common to the two proposed directions and is evaluated once. *)
initial40={run40["baseline",2,0.],run40["quadratic_minus",2,-.08],run40["quadratic_plus",2,.08],
 run40["quartic_minus",4,-.08],run40["quartic_plus",4,.08]};
pick40=Take[Ordering[initial40[[All,"accepted_cone_log_score"]],-2],2];
validated40=Table[Join[<|"selected_by_N64_score"->initial40[[k,"name"]],"N128_soft32"->run40[initial40[[k,"name"]],initial40[[k,"power"]],initial40[[k,"q"]],128,32]|>],{k,pick40}];
result40=<|"status"->"r040 synchronous curve-geometry deformation; floating direct split quadrature, not a certificate",
 "wolfram"->$Version,"n"->n,"base_parameters"-><|"h"->h40,"eta"->eta40,"source_scale"->scale40|>,
 "source_weight_formula"->"w_j exp(.2(v_j/16)^2) sum(w)/sum(w exp(.2(v_j/16)^2))",
 "geometry_formula"->"d_j -> d_j + q (v_j/16)^power, applied to both source shifts and receiver anchors",
 "initial_N64_soft32"->initial40,"selected_N128_soft32"->validated40,"scalar_reference_log_score"->71.71499208794934,
 "selection"->"two largest numerical accepted-cone scores at N64; all accepted groups require average > 1+1e-5"|>;
Export[out40,result40,"RawJSON"];Print[<|"N64"->initial40,"N128"->validated40|>];
