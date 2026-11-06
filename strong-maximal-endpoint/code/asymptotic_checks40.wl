(* Exact and high-precision checks for r040. Not interval certificates. *)
here40=DirectoryName[$InputFileName];args40=$ScriptCommandLine;
out40=If[Length[args40]>=2,args40[[2]],FileNameJoin[{$TemporaryDirectory,"asymptotic_checks40.json"}]];
$rate39LibraryMode=True;Get[FileNameJoin[{here40,"global_inflation_check39.wl"}]];$rate39LibraryMode=False;
A340=mu(mu-1)/8;K440=mu(mu-1)(3+5mu)/384;
th40[e_?NumericQ]:=z/.FindRoot[lsp[z]==yp[e],{z,e-e^2/10},WorkingPrecision->(wp39-10),AccuracyGoal->(wp39-35),PrecisionGoal->(wp39-35)];
row40[ep_,cc_]:=Module[{e=N[cc ep,wp39+20],th,y,ist,j,jp,ff,hh},y=yp[e];th=th40[e];ist=th y-ls[th];ff=ii[e]-ist-dd[e];j=ist-ep dd[e];jp=th-ep de[e]/ye[e];hh=j-y jp+ls[jp];<|"eps"->N[ep,18],"eta"->N[e,18],"F"->N[ff,24],"Jprime"->N[jp,24],"Hnew"->N[hh,24]|>];
rateRows40=Flatten[Table[row40[ep,cc],{ep,{1/100,1/200,1/500,1/1000}},{cc,{3,6}}],1];
crossC40=(1/(6(1/100)A340))^(1/3);
millsRows40=Table[Module[{eta=N[crossC40(Log[tt]/tt)^(1/3),wp39+20],mg},mg=Sqrt[Pi mu tt/2]Exp[mu tt eta^2/2]Erfc[eta Sqrt[mu tt/2]];<|"t"->tt,"eta"->N[eta,20],"scaled_Mills_error"->N[(eta mg-1)tt eta^2,24],"limit"->N[-1/mu,24]|>],{tt,{10^6,10^9,10^12}}];
syms40=<|"normal_Mills_expansion"->ToString[Integrate[Exp[-u]Normal[Series[Exp[-u^2 v/(2nu)],{v,0,3}]],{u,0,Infinity},Assumptions->nu>0],InputForm],
 "fixed_cutoff_remaining_exponential_mass_limit"->ToString[Limit[Exp[-aa Sqrt[Log[tt]/tt]cc Log[tt]],tt->Infinity,Assumptions->aa>0&&cc>0],InputForm],
 "Fourier_second_derivative_Gaussian_integral"->ToString[Integrate[u^2 Exp[-cc tt u^2],{u,-Infinity,Infinity},Assumptions->cc>0&&tt>0],InputForm],
 "optimized_Bernoulli_tail_factor"->ToString[Limit[(1+Exp[-Sqrt[d]])/(1+Exp[-1/Sqrt[d]]),d->0,Direction->"FromAbove"],InputForm],
 "fixed_k_minimal_tail_enhancements"->Table[<|"k"->kk,"enhancement_at_alpha_1"->N[(p/8)(HarmonicNumber[kk]-1)/(kk+1),24]|>,{kk,{4,8,16}}]|>;
Export[out40,<|"wolfram"->$Version,"base_working_precision"->wp39,"A3"->N[A340,24],"K4"->N[K440,24],"Ccrit_eps_0_01"->N[crossC40,24],"rate_rows"->rateRows40,"Gaussian_Mills_model_only"->millsRows40,"symbolic_checks"->syms40|>,"RawJSON"];
Print[<|"output"->out40,"rate_cases"->Length[rateRows40],"base_precision"->wp39|>];
