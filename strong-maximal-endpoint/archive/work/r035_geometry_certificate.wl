(* Exact rational certificate: a true cross-component halo box outside every
   individual unit-component strong-maximal halo. *)
root=DirectoryName[$InputFileName];out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{root,"r035_geometry_certificate.json"}]];
eBox={{-4,-7/2},{-4,-7/2}}; comps={{{0,1},{1,2}},{{1,2},{2,3}},{{2,3},{1,2}},{{3,4},{2,3}}};amp=16;mass=64;
(* For x,y in E, [x,4]x[y,3] contains all four cells; its largest area is at
   (-4,-4), hence the displayed lower bound is uniform on the whole box. *)
fullMin=mass/((4-(-4))(3-(-4)));
(* For a single [a,b]x[c,d] unit component and receiver left/below it, the
   exact noncentered strong maximal value is amp/((b-x)(d-y)); its E maximum
   is attained at the upper-right E corner. *)
singleMax=Function[q,amp/((q[[1,2]]-(-7/2))(q[[2,2]]-(-7/2)))]/@comps;
mat={{0,1,0,0},{0,0,1,0},{0,1,0,0},{0,0,1,0}};
checks=<|"full_average_minus_one"->(fullMin-1),"all_single_below_one"->And@@Thread[singleMax<1],"max_single"->Max[singleMax],"rank"->MatrixRank[mat]|>;
If[!(fullMin>1&&And@@Thread[singleMax<1]&&MatrixRank[mat]==2),Print["certificate check failed"];Exit[1]];
res=<|"status"->"exact rational continuous cross-component halo certificate","input_cells"->comps,"amplitude"->amp,
 "receiver_box"->eBox,"receiver_box_area"->1/4,"full_witness_rectangle"->"[x,4] x [y,3] for each (x,y) in receiver_box","full_mass"->mass,"uniform_full_average_lower"->fullMin,
 "single_component_exact_maxima_on_box"->singleMax,"single_component_maximum"->Max[singleMax],"strict_full_margin"->fullMin-1,"strict_single_gap"->1-Max[singleMax],
 "binary_shape_matrix"->mat,"binary_shape_rank"->MatrixRank[mat],"conclusion"->"receiver_box lies in the full input strong halo but in none of the four individual unit-component strong halos"|>;
Export[out,res,"RawJSON"];Print[checks];
