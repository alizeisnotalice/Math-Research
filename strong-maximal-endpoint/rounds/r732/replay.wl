base=DirectoryName[$InputFileName];data=Import[FileNameJoin[{base,"capacity.json"}],"RawJSON"];
atomize[boxes_]:=Module[{k=Length[boxes],sub,iv,av,keep,inc,mass},sub=Table[Select[Range[k],BitAnd[mask,2^(#-1)]>0&],{mask,1,2^k-1}];iv=Table[Times@@Table[Max[0,Min[boxes[[ss,j,2]]]-Max[boxes[[ss,j,1]]]],{j,Length[First[boxes]]}],{ss,sub}];av=iv;Do[av[[i]]=iv[[i]]-Sum[If[BitAnd[i,j]==i,av[[j]],0],{j,i+1,2^k-1}],{i,2^k-1,1,-1}];keep=Select[Range[Length[av]],av[[#]]>0&];mass=av[[keep]];inc=Table[Boole[BitAnd[mask,2^(i-1)]>0],{i,k},{mask,keep}];{mass,inc,inc.mass,keep}];
rat[s_]:=ToExpression[s];
dec[s_]:=Module[{ss=StringReplace[s,"-"->""],sg=If[StringStartsQ[s,"-"],-1,1],parts},parts=StringSplit[ss,"."];sg FromDigits[StringJoin[parts]]/10^StringLength[Last[parts]]];
check[r_]:=Module[{mass,inc,vol,masks,om,pi,nu,act,pir,am,pn,q,g,ff,lower,upper,hall,f,ell,av,price,density,remain,npprice=0,ids,k,primal,numeric},
 {mass,inc,vol,masks}=atomize[r["boxes"]];om=Total[mass];pi=rat/@r["pi"];nu=rat/@r["nu"];act=r["active_rows"]+1;pir=pi[[act]];k=Length[pi];
 hall=And@@Table[Total[Pick[pi,Table[Boole[BitAnd[s,2^(i-1)]>0],{i,k}],1]]<=Total[Pick[mass,Boole[BitAnd[s,#]>0]&/@masks,1]]/om,{s,1,2^k-1}];
 f=rat/@r["source_f"];ell=rat/@r["source_ell"];av=inc.(mass f)/vol;price=rat[r["exact_NP_price"]];
 primal=And@@Join[Thread[f>=0],Thread[f<=1],Thread[ell>=0],Thread[ell<=1],Flatten[Table[If[inc[[i,j]]==1,{ell[[j]]<=av[[i]]},{}],{i,k},{j,Length[mass]}]],{mass.ell>=rat[r["target"]]om,mass.f==price}];
 density=Transpose[inc].(pi/vol);remain=rat[r["target"]];ids=Reverse[Ordering[density]];
 Do[If[remain<=mass[[j]]density[[j]],npprice+=remain/density[[j]];remain=0;Break[],npprice+=mass[[j]];remain-=mass[[j]]density[[j]]],{j,ids}];
 am=Table[mass inc[[i]]/vol[[i]],{i,act}];pn=nu.am;q=mass/om;g=am.Log[pn/q]-Log[nu/pir];ff=nu.g;
 lower=dec[r["capacity_lower"]];upper=dec[r["capacity_upper"]];numeric=TrueQ[N[ff-lower,80]>=0&&N[upper-Max[g],80]>=0];
 <|"case"->r["case"],"exact_geometry"->(mass==rat/@r["atom_masses"]&&vol==rat/@r["row_volumes"]&&om==rat[r["union_volume"]]),"exact_Hall"->TrueQ[hall&&Total[pi]==1&&Min[pi]>=0],"exact_primal_dual_price"->TrueQ[primal&&remain==0&&npprice==price],"positive_candidate"->TrueQ[Min[nu]>0&&Total[nu]==1],"capacity_80digit_within_rational_certificate"->numeric,"F80"->ToString[N[ff,80],InputForm],"upper_gradient80"->ToString[N[Max[g],80],InputForm]|>
];
LaunchKernels[12];DistributeDefinitions[atomize,rat,dec,check];res=ParallelMap[check,data["results"]];CloseKernels[];
Export[FileNameJoin[{base,"replay_wolfram.json"}],<|"status"->"complete","kernels"->12,"cases"->Length[res],"scope"->"Independent exact geometry/Hall/source price replay and80digit capacity evaluation, not an independent interval-log algorithm.","results"->res|>,"RawJSON"];
Print[<|"cases"->Length[res],"all_pass"->And@@Flatten[Lookup[res,{"exact_geometry","exact_Hall","exact_primal_dual_price","positive_candidate","capacity_80digit_within_rational_certificate"}]]|>];Quit[];
