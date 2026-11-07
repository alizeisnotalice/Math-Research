(* Independent local check of a nonuniform outer mesh. Numeric corroboration only. *)
ClearAll[phlog,calc];
phlog[n_,x_]:=If[x<=1,Log[x],Log[x]+(n-1)Log[1+Log[x]]];
calc[n_Integer,prec_Integer]:=Module[{m=Floor[n^(1/4)],q=Ceiling[N[Log[n]^2,100]],r,al,w1,qq,rr,eta,bs,bsStrict,low1,low,high1,high,masses,logterms,logsum,raw,scaled,strict,pcost},
 r=q/n;al=2r;w1=1-(m-1)al;qq=1+n/2;rr=qq^(m-1);eta=r(1-w1)/(w1-r);bs=1+eta;bsStrict=bs(1+1/n^2);
 low1=1-r/w1;low=1/2-1/n;high1=1+r/(w1(qq-1));high=1+1/n;
 masses={w1/rr,al(1-1/rr),w1(qq-1)/rr,al(1-qq/rr)};
 pcost[z_]:=Module[{hh={low1,low,high1,high},ma=Append[masses,1-Total[masses]],xs,logs,mm},xs=Append[hh,1];logs=MapThread[If[#1==0,-Infinity,N[Log[#1]+phlog[n,z #2],prec]]&,{ma,xs}];mm=Max[logs];mm+Log[Total[Exp[logs-mm]]]];
 raw=pcost[1];scaled=pcost[bs];strict=pcost[bsStrict];
 <|"N"->n,"M"->m,"q"->q,"precision"->prec,"feasible"->TrueQ[w1>r&&low>0&&Min[masses]>=0],"eta"->N[eta,prec],"N_eta"->N[n eta,prec],"log_raw_cost"->N[raw,prec],"log_scaled_cost"->N[scaled,prec],"log_strict_scaled_cost"->N[strict,prec],"strong_at_strict_scale"->N[(1+1/n^2)(1-r),prec],"adaptive_at_strict_scale"->N[1+1/n^2,prec],"branch_equality"->TrueQ[(1-r)/(1+eta)==(1-w1)/(1-w1+eta)]|>];

small=Flatten[Table[Module[{r=1/(10 mm),al,w,qq,rs,hh,lev,prefix,eta,glob,claimed,recMass},
 al=2r;w=Join[{1-(mm-1)al},ConstantArray[al,mm-1]];qq=1+nn/2;
 rs=Join[{0},Table[qq^(j-1),{j,mm}]];hh=Differences[rs];
 lev=Table[1+If[i==j,-r rs[[i+1]]/(w[[i]]hh[[i]]),0]+If[i+1==j,r rs[[i+1]]/(w[[i]]hh[[i+1]]),0],{i,mm},{j,mm}];
 prefix=Table[Sum[hh[[k]]lev[[i,k]],{k,j}]/rs[[j+1]],{i,mm},{j,mm}];
 eta=r(1-w[[1]])/(w[[1]]-r);
 glob=Max[Flatten[Table[Sum[w[[i]]prefix[[i,j]],{i,s,mm}]/(Total[w[[s;;mm]]]+eta),{s,mm},{j,mm}]]];
 claimed=(1-r)/(1+eta);
 recMass=Sum[w[[i]]hh[[i+1]]/Last[rs],{i,mm-1}];
 <|"N"->nn,"M"->mm,"prefix_check"->TrueQ[prefix==Table[1-If[i==j,r/w[[i]],0],{i,mm},{j,mm}]],"aggregate_cell_check"->TrueQ[And@@Thread[w.lev==ConstantArray[1-r,mm]]],"global_suffix_check"->TrueQ[glob==claimed],"recovery_mass_check"->TrueQ[recMass==w[[1]](qq-1)/qq^(mm-1)+al(1-qq/qq^(mm-1))]|>],{nn,{4,16,64}},{mm,{2,3,5}}],1];
If[!And@@Flatten[Lookup[small,#]&/@{"prefix_check","aggregate_cell_check","global_suffix_check","recovery_mass_check"}],Exit[3]];

rows=calc[#,70]&/@{2^10,2^16,2^24,2^32,2^40};
rows2=calc[#,140]&/@{2^10,2^16,2^24,2^32,2^40};
diff=Table[<|"N"->rows[[i]]["N"],"difference"->ToString[InputForm[Abs[rows[[i]]["log_strict_scaled_cost"]-rows2[[i]]["log_strict_scaled_cost"]]]],"pass"->TrueQ[Abs[rows[[i]]["log_strict_scaled_cost"]-rows2[[i]]["log_strict_scaled_cost"]]<10^-60]|>,{i,Length[rows]}];
If[!And@@Lookup[rows,"feasible"]||!And@@Lookup[rows,"branch_equality"]||!And@@Lookup[diff,"pass"],Exit[2]];
out=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/weighted_check60.json"];
Export[out,<|"wolfram"->$Version,"small_exact"->small,"rows"->rows,"precision_doubling"->diff,"cube_halo_limit_upper"->N[BesselI[0,2 Sqrt[2]],70],"scope"->"Weighted pointwise mesh; not an endpoint lower bound"|>,"RawJSON"];
Print[KeyTake[#,{"N","M","N_eta","log_scaled_cost","log_strict_scaled_cost"}]&/@rows];
