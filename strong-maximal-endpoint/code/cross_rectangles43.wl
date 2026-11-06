(* r043: exact finite-grid cross-component rectangle scan in dimension 2.
   A counted receiving cell is wholly covered by an actual strict-average
   rectangle.  Thus every reported union is a genuine continuous witness
   subset, complete over the stated finite edge family (not the full halo). *)
root43=DirectoryName[$InputFileName];out43=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{root43,"r043_cross_rectangles.json"}]];
tau43=1001/1000;edges43={-8,-4,-1,0,1,2,3,4,5,8,12};ne43=Length[edges43]-1;
inBit43[i_,j_]:=2^(i-1+4(j-1));outBit43[i_,j_]:=2^(i-1+ne43(j-1));
pop43[x_Integer]:=DigitCount[x,2,1];
mask43[cells_]:=Total[(inBit43@@#)&/@cells];
area43[m_Integer]:=Total[Flatten@Table[If[BitAnd[m,outBit43[i,j]]>0,(edges43[[i+1]]-edges43[[i]])(edges43[[j+1]]-edges43[[j]]),0],{i,ne43},{j,ne43}]];
rect43=Flatten[Table[Module[{xl=edges43[[a]],xr=edges43[[b]],yl=edges43[[c]],yr=edges43[[d]],im=0,om=0,ar},ar=(xr-xl)(yr-yl);
 Do[If[xl<=i-1&&i<=xr&&yl<=j-1&&j<=yr,im+=inBit43[i,j]],{i,4},{j,4}];
 Do[If[xl<=edges43[[i]]&&edges43[[i+1]]<=xr&&yl<=edges43[[j]]&&edges43[[j+1]]<=yr,om+=outBit43[i,j]],{i,ne43},{j,ne43}];
 <|"input_mask"->im,"output_mask"->om,"area"->ar,"bounds"->{xl,xr,yl,yr}|>],{a,ne43},{b,a+1,ne43+1},{c,ne43},{d,c+1,ne43+1}],3];
halo43[m_Integer,a_Integer]:=Module[{good=Select[rect43,a pop43[BitAnd[m,#"input_mask"]]/#"area">tau43&]},
 <|"mask"->Fold[BitOr,0,good[[All,"output_mask"]]],"good_rectangles"->Length[good],"least_margin"->If[Length[good]>0,Min[a pop43[BitAnd[m,#"input_mask"]]/#"area"-tau43&/@good],-Infinity]|>];
phi243[a_Integer]:=a(1+Log[a]);
record43[name_String,parts_List,a_Integer]:=Module[{masks=mask43/@parts,m,all,each,comp,cmask,cross,cost},m=Total[masks];all=halo43[m,a];each=halo43[#,a]&/@masks;
 cmask=Fold[BitOr,0,each[[All,"mask"]]];cross=BitAnd[all["mask"],BitNot[cmask]];cost=pop43[m] phi243[a];
 <|"pattern"->name,"components"->parts,"amplitude"->a,"input_cell_count"->pop43[m],
 "actual_orlicz_denominator_formula"->ToString[pop43[m] a (1+Log[a]),InputForm],"actual_orlicz_denominator_numeric"->N[cost,20],
 "finite_family_rectangle_count"->Length[rect43],"full_witness_union_area"->area43[all["mask"]],"component_witness_union_area"->area43[cmask],
 "cross_component_witness_area"->area43[cross],"cross_area_fraction_of_full"->N[area43[cross]/area43[all["mask"]],16],
 "full_good_rectangle_count"->all["good_rectangles"],"component_good_rectangle_counts"->each[[All,"good_rectangles"]],
 "least_full_rectangle_margin_over_tau"->ToString[all["least_margin"],InputForm],
 "full_ratio"->N[area43[all["mask"]]/cost,20],"component_only_ratio"->N[area43[cmask]/cost,20],
 "cross_ratio_gain"->N[(area43[all["mask"]]-area43[cmask])/cost,20]|>];
patterns43=<|
 "horizontal_two_1x2"->{{{1,2},{1,3}},{{3,2},{3,3}}},
 "diagonal_two_1x2"->{{{1,1},{1,2}},{{3,3},{3,4}}},
 "diagonal_two_2x2"->{{{1,1},{1,2},{2,1},{2,2}},{{3,3},{3,4},{4,3},{4,4}}},
 "offset_two_1x3"->{{{1,1},{1,2},{1,3}},{{3,2},{3,3},{3,4}}}
 |>;amps43={2,3,4,8,16};rows43=Flatten[KeyValueMap[Table[record43[#1,#2,a],{a,amps43}]&,patterns43],1];
bestFull43=First@SortBy[rows43,-#"full_ratio"&];bestGain43=First@SortBy[rows43,-#"cross_ratio_gain"&];
res43=<|"status"->"r043 exact finite-grid full-vs-component rectangle unions; no continuum optimum or high-dimensional lift claimed",
 "dimension"->2,"threshold"->"1001/1000","input_grid"->"[0,4]^2 unit cells","receiver_edges"->edges43,
 "scope"->"For each row, union every valid grid-edge rectangle; each counted output cell is fully contained in a strict-average rectangle, hence is a true continuous witness. Completeness is only for these 3025 rectangles and these receiver cells.",
 "scan_count"->Length[rows43],"rows"->rows43,"best_full_ratio"->bestFull43,"largest_cross_ratio_gain"->bestGain43,
 "comparison"->"component_witness_union uses the same complete finite rectangle family applied separately to the two original boxes. It is not asserted to equal either component's full continuous halo.",
 "interpretation"->"Actual Phi_2 denominators are reported for each same input; cross_ratio_gain is geometric finite-family gain only and is not a C2 or scalar-constant comparison."|>;
Export[out43,res43,"RawJSON"];Print[<|"best_full"->KeyTake[bestFull43,{"pattern","amplitude","full_ratio","component_only_ratio","cross_ratio_gain"}],"best_cross"->KeyTake[bestGain43,{"pattern","amplitude","full_ratio","component_only_ratio","cross_ratio_gain"}]|>];
