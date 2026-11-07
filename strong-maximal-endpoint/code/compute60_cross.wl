(* Independent native Wolfram evaluation of the final SciPy-selected source.
   Fixed source weights; no optimization or multiplier calibration. *)
args60=$ScriptCommandLine;repo60=DirectoryName[DirectoryName[$InputFileName]];
in60=If[Length[args60]>=2,args60[[2]],FileNameJoin[{repo60,"rounds","r060","compute60.json"}]];
out60=If[Length[args60]>=3,args60[[3]],"/tmp/compute60_cross.json"];
Block[{$ScriptCommandLine={},$r032LibraryMode=True},Get[FileNameJoin[{repo60,"code","receiver_continuous32.wl"}]]];
s60=Import[in60,"RawJSON"];par60=s60["parameters"];row60=SelectFirst[s60["results"],#["gamma"]==s60["best_by_128node_score"]&];
n=par60["n"];b=par60["b"];cap=par60["cap"];c=N[par60["c"]];tau=par60["tau"];w0=row60["fixed_weight_refinements"]["corrected_weights"];
vv=N[par60["v_nonnegative"]];dd=row60["gamma"]vv^2/n;ng=Length[w0];coeff=Join[{1.},ConstantArray[2.,ng-1]];
u=Join[{{0.,0.}},Transpose[{Rest[dd]+Rest[vv],Rest[dd]-Rest[vv]}],Transpose[{Rest[dd]-Rest[vv],Rest[dd]+Rest[vv]}]];
lo=First[row60["finite_window"]];hi=Last[row60["finite_window"]];pc=buildPolyCV[b,c];
Clear[probVecC];probVecC[t1_?NumericQ,t2_?NumericQ,cc_List]:=Module[{full},full=polyVec[t1-u[[All,1]]]polyVec[t2-u[[All,2]]]Exp[Total[u,{2}]-t1-t2];Join[{First[full]},full[[2;;ng]]+full[[ng+1;;2ng-1]]]];
start60=AbsoluteTime[];res60=evaluateReceiverC[ConstantArray[c,ng],1.,128,32];ref60=Last[row60["fixed_weight_refinements"]["rows"]];
checks60=<|"logJ_difference"->res60["actual_logJ"]-ref60["logJ"],"all_average_max_difference"->Max[Abs[res60["averages_by_group"]-ref60["allaverages"]]]|>;
If[!TrueQ[Abs[checks60["logJ_difference"]]<10^-7&&checks60["all_average_max_difference"]<10^-7],Print[checks60];Exit[2]];
Export[out60,<|"status"->"fixed-source independent floating cross-check passed","wolfram"->$Version,"seconds"->N[AbsoluteTime[]-start60],"result"->res60,"checks"->checks60,"scope"->"Branch-split finite-window floating quadrature, not interval certification"|>,"RawJSON"];Print[checks60];
