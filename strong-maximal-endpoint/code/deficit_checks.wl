root=DirectoryName[DirectoryName[$InputFileName]];
Clear[aa,bb,pp,ss,cc];
primitive=FullSimplify[pp(pp-1)Integrate[ss^(pp-2)(ss-aa),{ss,aa,bb},
 Assumptions->{bb>aa>0,pp>1}]-((pp-1)bb^pp-pp aa bb^(pp-1)+aa^pp),
 Assumptions->{bb>aa>0,pp>1}];
quadratic=Table[FullSimplify[D[(p-1)cc^p-p cc^(p-1)+1,{cc,2}]-p(p-1),Assumptions->cc>=1],{p,{2,4,8}}];
strongCandidates={20/9,5/3,5/3,5/2};
first=Maximize[{(25/3-(10/3)cc)/(3-cc),0<=cc<=1},cc];
second=Maximize[{5(2-cc)/(3-cc),1<=cc<=2},cc];
composition=Max[first[[1]],second[[1]],0]; strong=Max[strongCandidates];
result=<|"wolfram_version"->$Version,"overlap_primitive_residual"->ToString[primitive,InputForm],
 "quadratic_second_derivative_slacks"->(ToString[#,InputForm]&/@quadratic),
 "strong_rectangle_at_3_3"->ToString[strong,InputForm],
 "coordinate_composition_at_3_3"->ToString[composition,InputForm],
 "strict_gap"->ToString[composition-strong,InputForm],
 "annulus_true_max_at_5_4"->ToString[Max[{4/9,2/7,2/3,0}],InputForm],
 "annulus_anchored_max_at_5_4"->ToString[Max[{4/9,2/7,2/5}],InputForm]|>;
Export[FileNameJoin[{root,"rounds","r010","symbolic_checks.json"}],result,"RawJSON"];
Print[ExportString[result,"RawJSON"]];
