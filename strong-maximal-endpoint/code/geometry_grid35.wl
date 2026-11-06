(* r035: exact finite-grid actual-rectangle search for disconnected 4-cell inputs.
   Every accepted receiver cell is contained in one explicitly enumerated rectangle
   whose exact average is > tau.  This is a finite-grid lower witness only. *)
root=DirectoryName[$InputFileName];out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{$TemporaryDirectory,"geometry_grid35.json"}]];mode=If[Length[$ScriptCommandLine]>=3,$ScriptCommandLine[[3]],"full"];
tauNum=1001;tauDen=1000;edges={-4,-1,0,1,2,3,4,5,8};ne=Length[edges]-1;
cellBit[i_,j_]:=2^(i+ne (j-1)); inBit[i_,j_]:=2^(i-1+4 (j-1));
pop[x_Integer]:=DigitCount[x,2,1];
outArea[mask_Integer]:=Total[Flatten[Table[If[BitAnd[mask,cellBit[i,j]]>0,(edges[[i+1]]-edges[[i]])(edges[[j+1]]-edges[[j]]),0],{i,1,ne},{j,1,ne}]]];
rects=Flatten[Table[Module[{xl=edges[[a]],xr=edges[[b]],yl=edges[[c]],yr=edges[[d]],im=0,om=0},
  Do[If[xl<=i-1&&i<=xr&&yl<=j-1&&j<=yr,im+=inBit[i,j]],{i,1,4},{j,1,4}];
  Do[If[xl<=edges[[i]]&&edges[[i+1]]<=xr&&yl<=edges[[j]]&&edges[[j+1]]<=yr,om+=cellBit[i,j]],{i,1,ne},{j,1,ne}];
  {im,om,(xr-xl)(yr-yl),{xl,xr,yl,yr}}],{a,1,ne},{b,a+1,ne+1},{c,1,ne},{d,c+1,ne+1}],3];

halo[mask_Integer,amp_Integer]:=Module[{u=0,cnt},Do[cnt=pop[BitAnd[mask,r[[1]]]];If[tauDen amp cnt>tauNum r[[3]],u=BitOr[u,r[[2]]]],{r,rects}];u];
components[mask_Integer]:=Module[{left=mask,ans={},seed,front,nb,idx,i,j},
 While[left>0,seed=First@FirstPosition[IntegerDigits[left,2,16],1]; (* index is MSB-order, so rebuild safely below *)
  seed=First@Select[Range[16],BitAnd[left,2^(#-1)]>0&];front={seed};nb=0;
  While[Length[front]>0,idx=First[front];front=Rest[front];If[BitAnd[nb,2^(idx-1)]==0,
    nb=BitOr[nb,2^(idx-1)];i=Mod[idx-1,4]+1;j=Quotient[idx-1,4]+1;
    front=Join[front,Select[{If[i>1,idx-1,0],If[i<4,idx+1,0],If[j>1,idx-4,0],If[j<4,idx+4,0]},#>0&&BitAnd[left,2^(#-1)]>0&]]]];
  AppendTo[ans,nb];left=BitAnd[left,BitNot[nb]]];ans];
coords[mask_Integer]:=Cases[Flatten[Table[If[BitAnd[mask,inBit[i,j]]>0,{{i,j}},{}],{i,1,4},{j,1,4}],2],{_Integer,_Integer}];
record[mask_Integer,amp_Integer]:=Module[{u,cs,uc,base,cost},u=halo[mask,amp];cs=components[mask];uc=Fold[BitOr,0,halo[#,amp]&/@cs];cost=N[pop[mask] amp (1+Log[amp]),30];
 <|"support"->coords[mask],"components"->(coords/@cs),"component_count"->Length[cs],"amplitude"->amp,"input_cells"->pop[mask],"actual_rectangle_union_area"->outArea[u],"finite_template_component_witness_union_area"->outArea[uc],"finite_template_cross_witness_increment"->outArea[u]-outArea[uc],"orlicz_cost"->cost,"grid_ratio"->N[outArea[u]/cost,20]|>];
subMasks=(Total[2^(#-1)]& /@ Subsets[Range[16],{4}]);If[mode==="quick",subMasks=Take[subMasks,20]];amps={2,4,8,16};
all=Flatten[Table[record[m,a],{m,subMasks},{a,amps}],1];disc=Select[all,#["component_count"]>=2&];
best=First@SortBy[disc,-#["grid_ratio"]&];bestCross=First@SortBy[Select[disc,#["finite_template_cross_witness_increment"]>0&],-#["finite_template_cross_witness_increment"]&];
uniform=Table[KeyTake[record[2^16-1,a],{"amplitude","input_cells","actual_rectangle_union_area","orlicz_cost","grid_ratio"}],{a,amps}];
res=<|"status"->"exact finite-grid actual-rectangle union search; finite family only, not C2 or high-dimensional certificate","tau"->"1001/1000","input_grid"->"[0,4]^2 divided into 4x4 unit cells","receiver_edges"->edges,"candidate_rectangles"->Length[rects],"mode"->mode,"support_count"->Length[subMasks],"search"->If[mode==="quick","first 20 four-cell supports; amplitudes 2,4,8,16","all 1820 four-cell supports; amplitudes 2,4,8,16"],"best_disconnected_ratio"->best,"largest_cross_union_increment"->bestCross,"uniform_4x4_templates"->uniform,
 "mechanism_note"->"union is formed over receiver cells covered by actual average-valid rectangles. The finite-template component witness union is only a diagnostic and is not the full true halo of a component."|>;
Export[out,res,"RawJSON"];Print[<|"best_ratio"->best["grid_ratio"],"best_cross_increment"->bestCross["finite_template_cross_witness_increment"],"best_support"->best["support"]|>];
