(* r047: exact finite-endpoint witnesses for the genuine (unrestricted)
   strong maximal halo of the r045 two-slab physical input.  This is a
   continuous lower subset of the full halo, not an equality computation. *)
root47=DirectoryName[$InputFileName];out47=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/fullhalo47.json"];
Clear[s47,e47,q47];beta47=1/4;f047[s_]=1-beta47(s-1/2)^2/2;g47[s_]=s-1/2;
fp47[ee_,qq_]:=f047[0]+ee qq g47[0];fm47[s_,ee_,qq_]:=f047[s]-beta47(1+s)(s-1/2)+ee qq(1+2s-1/2);
primM47[ee_,qq_,r_]:=Integrate[fm47[s47,ee,qq],{s47,0,r}];
over47[a_,b_,l_,r_]:={Max[a,l],Min[b,r]};
hmass47[xl_,xr_,ee_,qq_]:=Module[{pp,nn,ab},pp=over47[xl,xr,0,1];nn=over47[xl,xr,-1,0];
 If[pp[[2]]>pp[[1]],(pp[[2]]-pp[[1]])fp47[ee,qq],0]+If[nn[[2]]>nn[[1]],ab={-nn[[2]],-nn[[1]]};primM47[ee,qq,ab[[2]]]-primM47[ee,qq,ab[[1]]],0]];
massRect47[xl_,xr_,yl_,yr_,ee_]:=Sum[Module[{v=over47[yl,yr,sl[[1]],sl[[2]]]},If[v[[2]]>v[[1]],(v[[2]]-v[[1]])hmass47[xl,xr,ee,sl[[3]]],0]],{sl,{{0,1/2,1},{1/2,1,-1}}}];

(* Receiver cells and all rectangles use finite rational endpoint lists.
   A valid rectangle covers a cell only when it contains the entire cell. *)
xe47={-3,-2,-1,-1/2,0,1/2,1,3/2,2,3};ye47={-2,-1,0,1/4,1/2,3/4,1,2,3};nx47=Length[xe47]-1;ny47=Length[ye47]-1;
cells47=Flatten[Table[<|"ij"->{i,j},"bounds"->{xe47[[i]],xe47[[i+1]],ye47[[j]],ye47[[j+1]]},"area"->(xe47[[i+1]]-xe47[[i]])(ye47[[j+1]]-ye47[[j]])|>,{i,nx47},{j,ny47}],1];
rect47=Flatten[Table[Module[{xl=xe47[[a]],xr=xe47[[b]],yl=ye47[[c]],yr=ye47[[d]],inside},inside=Select[cells47,Function[cell,With[{bb=cell["bounds"]},xl<=bb[[1]]&&bb[[2]]<=xr&&yl<=bb[[3]]&&bb[[4]]<=yr]]];<|"bounds"->{xl,xr,yl,yr},"area"->(xr-xl)(yr-yl),"covered"->inside|>],{a,nx47},{b,a+1,nx47+1},{c,ny47},{d,c+1,ny47+1}],3];

cost47[ee_,prec_:90]:=Module[{extra=0,poly,roots,ed,lo,hi,mid,z},
 Do[poly=fm47[s47,ee,qq];roots=N[s47/.Solve[poly==1,s47,Reals],prec+20];ed=Join[{0},Sort@Select[roots,0<#<1&],{1}];
  Do[{lo,hi}=ed[[j;;j+1]];mid=(lo+hi)/2;If[(poly/.s47->mid)>1,extra+=NIntegrate[Evaluate[poly Log[poly]],{s47,lo,hi},WorkingPrecision->prec,AccuracyGoal->Floor[3prec/5],PrecisionGoal->Floor[prec/2],MaxRecursion->16]/2],{j,Length[ed]-1}];
  z=fp47[ee,qq];If[z>1,extra+=N[z Log[z],prec]/2],{qq,{1,-1}}];N[31/16+extra,prec]];
finiteHalo47[ee_]:=Module[{good,uc,ar,mm,one,covered},good=Select[rect47,Function[rr,massRect47[Sequence@@rr["bounds"],ee]/rr["area"]>1]];covered=If[Length[good]>0,Apply[Join,good[[All,"covered"]]],{}];uc=DeleteDuplicatesBy[covered,Function[cell,cell["ij"]]];ar=If[Length[uc]>0,Total[uc[[All,"area"]]],0];mm=If[Length[good]>0,Min[(massRect47[Sequence@@#["bounds"],ee]/#["area"]-1)&/@good],0];one=If[Length[good]>0,First[good],Missing["none"]];<|"valid_rectangle_count"->Length[good],"witness_cells"->Length[uc],"witness_cell_bounds"->If[Length[uc]>0,uc[[All,"bounds"]]],"witness_union_area_exact"->ToString[ar,InputForm],"least_rectangle_margin_exact"->ToString[mm,InputForm],"example_valid_rectangle"->If[one===Missing["none"],Null,one["bounds"]]|>];
row47[ee_]:=Module[{h=finiteHalo47[ee],cc=cost47[ee,100],area},area=ToExpression[h["witness_union_area_exact"]];Join[<|"epsilon"->ToString[ee,InputForm],"actual_Phi2_budget"->N[cc,30],"input_mass_exact"->"31/16"|>,h,<|"finite_witness_area_over_actual_budget"->N[area/cc,30]|>]];
eps47={0,1/32,1/16,1/8};rows47=row47/@eps47;
uniform47=<|"definition"->"(31/32) 1_{[-1,1]x[0,1]}","same_input_mass_exact"->"31/16","actual_Phi2_budget_exact"->"31/16","strict_full_halo_empty_reason"->"f<=31/32<1 pointwise, hence M_S f<1 everywhere","finite_witness_area"->0,"finite_witness_area_over_budget"->0|>;
res47=<|"status"->"r047 exact finite-endpoint full-rectangle witness unions for r045 input; not a full-halo-area computation",
 "input"->"r045 physical two-y-slab nonproduct f_epsilon supported on [-1,1]x[0,1]","threshold"->">1",
 "receiver_endpoint_grid"-><|"x"->xe47,"y"->ye47,"cells"->Length[cells47],"rectangles"->Length[rect47]|>,
 "method"->"All finite-endpoint rectangles are tested with exact rational polynomial integrals. A counted receiver cell lies wholly in a strict-average rectangle, so the union is a genuine subset of {M_S f_epsilon>1}. Phi2 budget is independently high-precision split integration at exact f=1 roots.",
 "rows"->rows47,"epsilon0_same_mass_reference"->First[rows47],"uniform_same_mass_baseline"->uniform47,
 "scope"->"The finite endpoint/cell family cannot equal the full continuous halo: arbitrary horizontal endpoints remain unoptimized. No numerical scan is presented as a certificate for full-halo area."|>;
Export[out47,res47,"RawJSON"];Print[rows47];
