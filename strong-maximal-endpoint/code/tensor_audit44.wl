(* r044: exact audit of tensor products of the fixed 20 r043 witnesses.
   This concerns only product receiving witnesses, never a full strong halo. *)
root44=DirectoryName[$InputFileName];in44=If[Length[$ScriptCommandLine]>=3,$ScriptCommandLine[[3]],FileNameJoin[{DirectoryName[root44],"rounds","r043","cross_rectangles43.json"}]];out44=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],FileNameJoin[{$TemporaryDirectory,"tensor_audit44.json"}]];
dat44=Import[in44,"RawJSON"];rows44=dat44["rows"];
trip44=({#["amplitude"],#["full_witness_union_area"],#["input_cell_count"],#["pattern"]}& /@ rows44);
ratio44[t_]:=t[[2]]/(t[[1]] t[[3]]);
amps44=Sort@DeleteDuplicates[trip44[[All,1]]];
bestAtAmp44=Association@Table[ToString[a]->First@SortBy[Select[trip44,#[[1]]==a&],-ratio44[#]&],{a,amps44}];
mAtAmp44=Association@Table[ToString[a]->ratio44[bestAtAmp44[ToString[a]]],{a,amps44}];
m44=Max[Values[mAtAmp44]]; If[m44!=4,Print["unexpected M"];Exit[1]];

(* Rigorous rational check: exp(69/100)<2 by a Taylor upper tail. *)
x44=69/100;nTaylor44=12;sTaylor44=Sum[x44^j/j!,{j,0,nTaylor44}];
rTaylor44=x44^(nTaylor44+1)/(nTaylor44+1)!/(1-x44/(nTaylor44+2)); eUpper44=sTaylor44+rTaylor44;
If[!(eUpper44<2),Print["log2 rational lower check failed"];Exit[1]];
ub44[k_Integer]:=m44^k/(1+k x44)^(2k-1);
ubRows44=Table[<|"k"->k,"upper_bound"->N[ub44[k],40],"log_upper_bound"->N[Log[ub44[k]],40]|>,{k,1,128}];
firstUB44=First[Select[ubRows44,#"upper_bound"<=1&]]["k"];

(* Finite-dictionary DP with 70-digit objective comparisons.  At a fixed amplitude only the greatest m_A
   can occur in an optimum.  A=2,4,8,16 contribute respectively 1,2,3,4
   copies of log 2; A=3 contributes one log 3. *)
m2=mAtAmp44["2"];m3=mAtAmp44["3"];m4=mAtAmp44["4"];m8=mAtAmp44["8"];m16=mAtAmp44["16"];
ww44={1,2,3,4};mm44={m2,m4,m8,m16};lm44=N[Log[mm44],70];lm3=N[Log[m3],70];l2=N[Log[2],70];l3=N[Log[3],70];
f44=ConstantArray[-Infinity,{129,513}];f44[[1,1]]=SetPrecision[0,70];
Do[Do[f44[[q+1,r+1]]=Max@Table[If[r>=ww44[[z]],f44[[q,r-ww44[[z]]+1]]+lm44[[z]],-Infinity],{z,4}],{r,0,4q}],{q,1,128}];
dp44[k_Integer]:=Module[{q,s,r,val,best={-Infinity,0,0,0}},
 Do[q=k-s;Do[If[f44[[q+1,r+1]]=!=-Infinity,
   val=f44[[q+1,r+1]]+s lm3-(2k-1) Log[1+r l2+s l3];
   If[val>best[[1]],best={val,q,r,s}]],{r,0,4q}],{s,0,k}];
 <|"log_product_witness"->N[best[[1]],40],"product_witness"->N[Exp[best[[1]]],40],"non3_factor_count"->best[[2]],"log2_weight"->best[[3]],"A3_factor_count"->best[[4]]|>
 ];
same44[k_Integer]:=Module[{a,val,best={-Infinity,0}},Do[val=Log[mAtAmp44[ToString[a]]^k]-(2k-1)Log[1+k Log[a]];If[val>best[[1]],best={val,a}],{a,amps44}];<|"log_same_witness"->N[best[[1]],40],"amplitude"->best[[2]],"same_witness"->N[Exp[best[[1]]],40]|>];
dpRows44=Table[Join[<|"k"->k|>,dp44[k],same44[k]],{k,1,128}];firstDP44=First[Select[dpRows44,#"product_witness"<=1&]]["k"];

(* Uniform proof: B_3<1 and B_(k+1)/B_k < 4/(1+(69/100)(k+1))^2<1 for k>=3. *)
b3=ub44[3];ratioBoundAt3=4/(1+x44*4)^2;
res44=<|"status"->"r044 exact tensor-product denominator audit over the fixed r043 finite witness dictionary; product receivers only",
 "wolfram"->$Version,"source_basename"->FileNameTake[in44],"source_sha256"->IntegerString[FileHash[in44,"SHA256"],16,64],"DP_comparison_precision"->70,"witness_count"->Length[rows44],
 "per_amplitude_best"->Association@KeyValueMap[(#1->With[{t=#2},<|"amplitude"->t[[1]],"pattern"->t[[4]],"V"->t[[2]],"support_area"->t[[3]],"m=V/(area*A)"->ToString[ratio44[t],InputForm]|>])&,bestAtAmp44],
 "M_exact"->ToString[m44,InputForm],"M_attained_by"->KeySelect[bestAtAmp44,mAtAmp44[#]==m44&],
 "product_formula"->"L=Product_i(V_i/(|E_i| A_i))/(1+Sum_i Log[A_i])^(2k-1)",
 "heterogeneous_upper_bound"->"L_(2k) <= 4^k/(1+k Log[2])^(2k-1) <= 4^k/(1+69 k/100)^(2k-1)",
 "log2_rational_lower_certificate"-><|"x"->"69/100","Taylor_order"->nTaylor44,"exp_upper_exact"->ToString[eUpper44,InputForm],"assertion"->"exp(69/100)<2, hence Log[2]>69/100"|>,
 "all_k_1_to_128_universal_bounds"->ubRows44,"first_universal_bound_le_1"->firstUB44,
 "uniform_decay_proof"-><|"B3_exact"->ToString[b3,InputForm],"B3_less_than_1"->(b3<1),"successive_ratio_upper_for_k_ge_3"->ToString[ratioBoundAt3,InputForm],
 "reason"->"B_(k+1)/B_k < 4/(1+(69/100)(k+1))^2, whose k=3 value is <1; therefore B_k<=B3<1 for every k>=3."|>,
 "log_decay"->"log L_(2k)/(2k) <= -log k + log(200/69)+1/2, hence -log k+O(1)",
 "finite_dictionary_DP_k_1_to_128"->dpRows44,"first_DP_product_witness_le_1"->firstDP44,
 "scope"->"No conclusion about full strong-maximal halos, arbitrary inputs, or tensor lifts beyond these explicitly product receiving witnesses."|>;
Export[out44,res44,"RawJSON"];Print[<|"M"->m44,"firstUB"->firstUB44,"firstDP"->firstDP44,"k1"->First[dpRows44],"k2"->dpRows44[[2]],"k128"->Last[dpRows44]|>];
