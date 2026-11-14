folder=DirectoryName[$InputFileName];
atomize[boxes_]:=Module[{k=Length[boxes],sub,iv,av,keep,inc,vol,mass},
 sub=Table[Select[Range[k],BitAnd[mask,2^(#-1)]>0&],{mask,1,2^k-1}];
 iv=Table[Times@@Table[Max[0,Min[boxes[[ss,j,2]]]-Max[boxes[[ss,j,1]]]],{j,Length[First[boxes]]}],{ss,sub}];av=iv;
 Do[av[[i]]=iv[[i]]-Sum[If[BitAnd[i,j]==i,av[[j]],0],{j,i+1,2^k-1}],{i,2^k-1,1,-1}];
 keep=Select[Range[Length[av]],av[[#]]>0&];mass=av[[keep]];inc=Table[Boole[BitAnd[mask,2^(i-1)]>0],{i,k},{mask,keep}];vol=inc.mass;{mass,inc,vol,keep}];
verify[boxes_,p_]:=Module[{mass,inc,vol,masks,k,x,ell,avg,cons,source,lam,z,t,sub,hall,price,true,om},
 {mass,inc,vol,masks}=atomize[boxes];k=Length[boxes];om=Total[mass];x=Array[xx,Length[mass]];ell=Array[ee,Length[mass]];
 avg=Table[Sum[inc[[i,j]]mass[[j]]x[[j]],{j,Length[mass]}]/vol[[i]],{i,k}];
 cons=Flatten[Table[If[inc[[i,j]]==1,{ell[[j]]<=avg[[i]]},{}],{i,k},{j,Length[mass]}]];
 source=First[Minimize[{mass.x,And@@Join[Thread[x>=0],Thread[x<=1],Thread[ell>=0],Thread[ell<=1],cons,{mass.ell>=p om}]},Join[x,ell]]];
 lam=Array[ll,k];z=Array[zz,Length[mass]];t=vol.lam;
 hall=Table[Total[Table[If[BitAnd[mask,2^(i-1)]>0,vol[[i]]lam[[i]],0],{i,k}]]<=Total[Table[If[BitAnd[mask,masks[[j]]]>0,mass[[j]],0],{j,Length[mass]}]]/om*t,{mask,1,2^k-1}];
 price=First[Maximize[{p t-mass.z,And@@Join[Thread[lam>=0],Thread[z>=0],Thread[z>=Transpose[inc].lam-1],hall]},Join[lam,z]]];
 true=First[Minimize[{mass.x,And@@Join[Thread[x>=0],Thread[x<=1],Thread[avg>=p]]},x]];
 <|"dimension"->Length[First[boxes]],"boxes"->boxes,"rows"->k,"atoms"->Length[mass],"p"->ToString[p,InputForm],"source_envelope"->ToString[source,InputForm],"responsible_Hall_dual"->ToString[price,InputForm],"true_source"->ToString[true,InputForm],"identity"->(source==price),"price_order"->(source<=true)|>];
recs=Select[Import[FileNameJoin[{folder,"..","r710","catalogues.json"}],"RawJSON"],#["n"]<=3&&#["K"]==2&];
boxes=Join[(Take[#["boxes"],Min[6,Length[#["boxes"]]]]&/@recs),{{{{0,10}},{{-1,3}},{{7,11}}}}];
cases=Flatten[Table[{b,p},{b,boxes},{p,{1/4,1/2,3/4}}],1];LaunchKernels[12];DistributeDefinitions[atomize,verify];res=ParallelMap[Function[c,TimeConstrained[MemoryConstrained[verify@@c,268435456,$Failed],120,$Aborted]],cases,Method->"CoarsestGrained"];CloseKernels[];
Export[FileNameJoin[{folder,"wolfram.json"}],<|"version"->$Version,"results"->res|>,"RawJSON"];Print[<|"cases"->Length[res],"identity_passes"->Count[res,r_Association/;TrueQ[r["identity"]&&r["price_order"]]]|>];Quit[];
