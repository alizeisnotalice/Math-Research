base=DirectoryName[$InputFileName];
source=Import[FileNameJoin[{base,"..","r732","capacity.json"}],"RawJSON"]["results"];
hulls=Association[(#["case"]->#)&/@Import[FileNameJoin[{base,"posterior_hull_capacity.json"}],"RawJSON"]["results"]];
gains=Association[(#["case"]->#)&/@Import[FileNameJoin[{base,"tilt_gain.json"}],"RawJSON"]["results"]];
rat[s_]:=ToExpression[s];
dec[s_]:=Module[{ss=StringReplace[s,"-"->""],sg=If[StringStartsQ[s,"-"],-1,1],p},p=StringSplit[ss,"."];sg FromDigits[StringJoin[p]]/10^StringLength[Last[p]]];
kl[p_,q_]:=Total[MapThread[If[#1==0,0,#1 Log[#1/#2]]&,{p,q}]];
elog[x_]:=elog[x]=Total[(Last[#]Log[First[#]]&)/@FactorInteger[x]];
expandlog[x_]:=Expand[x/.{Log[t_Rational]:>elog[t],Log[t_Integer]:>elog[t]}];
check[r_]:=Module[{mass,vol,masks,pi,act,ker,q,pm,alpha,k,a,pa,mi,res,ident,h,gn,z,nu,pn,g,ff,hh,hcheck,step,phi,gainvalue,oldg,analytic,gcheck},
 mass=rat/@r["atom_masses"];vol=rat/@r["row_volumes"];masks=r["atom_masks"];act=r["active_rows"]+1;pi=(rat/@r["pi"])[[act]];q=mass/Total[mass];k=Length[act];a=Length[mass];
 ker=Table[Table[If[BitAnd[masks[[j]],2^(i-1)]>0,mass[[j]]/vol[[i]],0],{j,a}],{i,act}];pm=pi.ker;alpha=Table[pi ker[[All,j]]/pm[[j]],{j,a}];
 ident=And@@Table[pa=alpha[[z]].ker;mi=Sum[If[alpha[[z,i]]ker[[i,j]]==0,0,alpha[[z,i]]ker[[i,j]]Log[ker[[i,j]]/pa[[j]]]],{i,k},{j,a}];res=kl[pa,q]-kl[alpha[[z]],pi]-Log[pm[[z]]/q[[z]]]+mi;TrueQ[expandlog[res]==0],{z,a}];
 h=hulls[r["case"]];z=rat/@h["space_weights"];nu=z.alpha;pn=nu.ker;g=ker.Log[pn/q]-Log[nu/pi];ff=nu.g;hh=alpha.g;
 hcheck=TrueQ[nu==rat/@h["induced_nu"]&&N[ff-dec[h["hull_capacity_lower"]],80]>=0&&N[dec[h["hull_capacity_upper"]]-Max[hh],80]>=0&&N[dec[r["capacity_lower"]]-Max[hh]-dec[h["full_minus_hull_lower"]],80]>=0];
 gn=gains[r["case"]];z=rat/@gn["space_weights"];nu=z.alpha;pn=nu.ker;gainvalue=kl[pn,q]-kl[nu,pi]-kl[pm,q];step=rat[gn["step_t"]];phi=rat/@gn["phi"];oldg=ker.Log[pm/q];analytic=(nu-pi).oldg-step^2 Total[pm phi^2];
 gcheck=TrueQ[nu==rat/@gn["induced_nu"]&&z==pm(1+step phi)&&N[gainvalue-dec[gn["direct_gain_lower"]],80]>=0&&N[analytic-dec[gn["analytic_direction_gain_lower"]],80]>=0&&N[gainvalue-analytic,80]>=0];
 <|"case"->r["case"],"all_atom_posterior_identity_exact"->ident,"hull_positive_candidate_and80digit_bounds"->hcheck,"explicit_gain_and80digit_bound"->gcheck,"explicit_gain80"->ToString[N[gainvalue,80],InputForm],"hull_gap_lower_check80"->ToString[N[dec[r["capacity_lower"]]-Max[hh],80],InputForm]|>
];
LaunchKernels[12];DistributeDefinitions[hulls,gains,rat,dec,kl,elog,expandlog,check];res=ParallelMap[check,source];CloseKernels[];
Export[FileNameJoin[{base,"independent_wolfram.json"}],<|"status"->"complete","version"->$Version,"kernels"->12,"precision_digits"->80,"scope"->"Exact posterior identities, exact induced-label maps, independent80digit bound checks; rational strict inequalities are independently replayed by Python.","results"->res|>,"RawJSON"];
Print[<|"cases"->Length[res],"all_pass"->And@@Flatten[Lookup[res,{"all_atom_posterior_identity_exact","hull_positive_candidate_and80digit_bounds","explicit_gain_and80digit_bound"}]]|>];Quit[];
