(* Exact witnesses for the best r035 finite-grid cross-component candidate. *)
root=DirectoryName[$InputFileName];out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{root,"r035_geometry_witness.json"}]];
e={-4,-1,0,1,2,3,4,5,8};ne=Length[e]-1;ib[i_,j_]:=2^(i-1+4(j-1));ob[i_,j_]:=2^(i+ne(j-1));pop[x_Integer]:=DigitCount[x,2,1];
rects=Flatten[Table[Module[{xl=e[[a]],xr=e[[b]],yl=e[[c]],yr=e[[d]],im=0,om=0},Do[If[xl<=i-1&&i<=xr&&yl<=j-1&&j<=yr,im+=ib[i,j]],{i,1,4},{j,1,4}];Do[If[xl<=e[[i]]&&e[[i+1]]<=xr&&yl<=e[[j]]&&e[[j+1]]<=yr,om+=ob[i,j]],{i,1,ne},{j,1,ne}];{im,om,(xr-xl)(yr-yl),{xl,xr,yl,yr}}],{a,1,ne},{b,a+1,ne+1},{c,1,ne},{d,c+1,ne+1}],3];
amp=4;whole=ib[1,2]+ib[1,3]+ib[3,2]+ib[3,3];parts={ib[1,2]+ib[1,3],ib[3,2]+ib[3,3]};valid[mask_,r_]:=1000 amp pop[BitAnd[mask,r[[1]]]]>1001 r[[3]];
hu[mask_]:=Fold[BitOr,0,Select[rects,valid[mask,#]&][[All,2]]];full=hu[whole];single=Fold[BitOr,0,hu/@parts];
cross=Select[Flatten[Table[{i,j},{i,1,ne},{j,1,ne}],1],BitAnd[full,ob[#[[1]],#[[2]]]]>0&&BitAnd[single,ob[#[[1]],#[[2]]]]==0&];
witness[cell_]:=Module[{rr},rr=First@Select[rects,BitAnd[#[[2]],ob[cell[[1]],cell[[2]]]]>0&&valid[whole,#]&&Count[parts,p_/;BitAnd[p,#[[1]]]>0]>=2&];<|"receiver_cell"->{e[[cell[[1]]]],e[[cell[[1]]+1]],e[[cell[[2]]]],e[[cell[[2]]+1]]},"rectangle"->rr[[4]],"high_cells_captured"->pop[BitAnd[whole,rr[[1]]]],"exact_average"->(amp pop[BitAnd[whole,rr[[1]]]]/rr[[3]])|>];
res=<|"status"->"exact integer finite-grid cross-rectangle witnesses","tau"->"1001/1000","support"->{{1,2},{1,3},{3,2},{3,3}},"amplitude"->amp,"full_union_area"->48,"component_only_union_area"->38,"cross_increment_area"->10,"cross_receiver_cells"->(witness/@cross)|>;
Export[out,res,"RawJSON"];Print[<|"cross_cells"->Length[cross],"cross_area"->10|>];
