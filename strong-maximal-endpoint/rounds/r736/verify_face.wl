base=DirectoryName[$InputFileName];enum=Import[FileNameJoin[{base,"lf_vertices.json"}],"RawJSON"];
paths={FileNameJoin[{base,"..","r732","capacity.json"}]};
sourcePath=SelectFirst[paths,FileExistsQ[#]&&IntegerString[FileHash[#,"SHA256"],16,64]==enum["source_sha256"]&];
source=Import[sourcePath,"RawJSON"]["results"];lf=Association[(#["case"]->#)&/@enum["results"]];
atomize[boxes_]:=Module[{k=Length[boxes],sub,iv,av,keep,inc,mass},sub=Table[Select[Range[k],BitAnd[mask,2^(#-1)]>0&],{mask,1,2^k-1}];iv=Table[Times@@Table[Max[0,Min[boxes[[ss,j,2]]]-Max[boxes[[ss,j,1]]]],{j,Length[First[boxes]]}],{ss,sub}];av=iv;Do[av[[i]]=iv[[i]]-Sum[If[BitAnd[i,j]==i,av[[j]],0],{j,i+1,2^k-1}],{i,2^k-1,1,-1}];keep=Select[Range[Length[av]],av[[#]]>0&];mass=av[[keep]];inc=Table[Boole[BitAnd[mask,2^(i-1)]>0],{i,k},{mask,keep}];{mass,inc,inc.mass,keep}];
rat[s_]:=ToExpression[s];normal[row_]:=Module[{ints=Together[row](LCM@@Denominator[Together[row]])},ints/(GCD@@Abs[ints])];
check[r_]:=Module[{mass,inc,vol,masks,om,k,a,f,av,pp,eta,dens,q,cons,rec,verts,m,diff,rk,ann,rr,coords,proj,facets={},ss,pts,ns,c,b,vv,obj={},sols={},one,member,affineCount,all,started=AbsoluteTime[]},
 {mass,inc,vol,masks}=atomize[r["boxes"]];om=Total[mass];k=Length[vol];a=Length[mass];f=rat/@r["source_f"];av=inc.(mass f)/vol;q=rat[r["target"]];pp=Table[Unique["p"],{k}];eta=Unique["eta"];dens=Transpose[inc].(pp/vol);
 cons=Join[Thread[pp>=0],{Total[pp]==1,pp.av==q,0<=eta<=Max[1/vol]},Table[Total[Pick[pp,Table[Boole[BitAnd[s,2^(i-1)]>0],{i,k}],1]]<=Total[Pick[mass,Boole[BitAnd[s,#]>0]&/@masks,1]]/om,{s,1,2^k-1}],Table[Which[f[[j]]==0,dens[[j]]<=eta,f[[j]]==1,dens[[j]]>=eta,True,dens[[j]]==eta],{j,a}]];
 rec=lf[r["case"]];verts=(rat/@#)&/@rec["unique_pi"];m=Length[verts];
 member=And@@Table[TrueQ[And@@cons/.Join[Thread[pp->(rat/@v["pi"])],{eta->rat[v["threshold"]]}]],{v,rec["vertices"]}];
 If[m==1,rk=0;ann=IdentityMatrix[k],diff=Rest[verts]-ConstantArray[First[verts],m-1];rk=MatrixRank[diff];ann=NullSpace[diff]];
 Do[obj=Join[obj,{v.(pp-First[verts]),-v.(pp-First[verts])}],{v,ann}];affineCount=Length[obj];
 If[rk>0,
  rr=Select[RowReduce[diff],#!=ConstantArray[0,k]&];coords=Table[First[Select[Range[k],row[[#]]!=0&]],{row,rr}];proj=verts[[All,coords]];
  Do[pts=proj[[ss]];ns=If[rk==1,{{1}},NullSpace[Rest[pts]-ConstantArray[First[pts],rk-1]]];
   If[Length[ns]==1,c=First[ns];b=c.First[pts];vv=(c.#-b)&/@proj;
    If[Max[vv]<=0&&Min[vv]<0,AppendTo[facets,normal[Join[c,{b}]]],If[Min[vv]>=0&&Max[vv]>0,AppendTo[facets,normal[-Join[c,{b}]]]]]],
   {ss,Subsets[Range[m],{rk}]}];facets=DeleteDuplicates[facets];
  obj=Join[obj,Table[Most[fac].pp[[coords]]-Last[fac],{fac,facets}]];
 ];
 one[objective_]:=Module[{s=TimeConstrained[MemoryConstrained[Maximize[{objective,And@@cons},Join[pp,{eta}]],268435456,$Failed],30,$Aborted]},If[ListQ[s]&&Length[s]==2,First[s],s]];
 sols=one/@obj;all=member&&And@@(TrueQ[#==0]&/@sols);
 <|"case"->r["case"],"unique_pi_count"->m,"projected_affine_dimension"->rk,"listed_vertices_exact_LF"->member,"all_LF_pi_in_listed_convex_hull"->all,"affine_hull_tests"->affineCount,"complete_hull_facets"->Length[facets],"exact_linear_maxima"->(ToString[#,InputForm]&/@sols),"wall_seconds"->(AbsoluteTime[]-started)|>
];
LaunchKernels[12];DistributeDefinitions[lf,atomize,rat,normal,check];res=ParallelMap[check,source];CloseKernels[];
Export[FileNameJoin[{base,"verify_face_wolfram.json"}],<|"status"->"complete","version"->$Version,"kernels"->12,"source_sha256"->enum["source_sha256"],"method"->"Exact linear optimization proves every LFpoint lies in listed convex hull: verify affine hull and all facets, exhaustively generated from all affine-independent vertex subsets; listed vertices independently LF.","results"->res|>,"RawJSON"];
Print[<|"cases"->Length[res],"all_pass"->And@@Flatten[Lookup[res,{"listed_vertices_exact_LF","all_LF_pi_in_listed_convex_hull"}]]|>];Quit[];
