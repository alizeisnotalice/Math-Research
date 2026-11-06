(* r046: stable tilted-scale block-density kernel.  Floating numerical
   diagnostics only; no new maximal-operator candidate is asserted. *)
root46=DirectoryName[$InputFileName];out46=If[Length[$ScriptCommandLine]>=2,$ScriptCommandLine[[2]],"/tmp/block_kernel46.json"];
constWP46=360;mu46=N[1+ProductLog[1/E],constWP46];
p46[d_?NumericQ]:=1/(1+Exp[d]);u46[d_?NumericQ]:=d p46[d];m46[d_?NumericQ]:=1+u46[d];nu46[d_?NumericQ]:=1+d^2 p46[d](1-p46[d]);
ratio46[d_?NumericQ]:=m46[d]/d;
solveD46[b_Integer,c_?NumericQ,t_?NumericQ,wp_Integer]:=Module[{target=N[t/(b c),wp],lo,hi,mid,its,tol,d,rr,br,rel},
 (* Since 1<=m(d)<=2, the unique d solving m(d)/d=target lies in
    [1/target,2/target].  This remains valid for extreme target values. *)
 lo=1/target;hi=2/target;its=Ceiling[(wp-10) Log[10]/Log[2]];tol=N[10^(-(wp-20)),wp];
 Do[mid=(lo+hi)/2;If[ratio46[mid]>target,lo=mid,hi=mid],{its}];d=(lo+hi)/2;
 rr=Abs[t/c-b ratio46[d]];br=Max[Abs[t/c-b ratio46[lo]],Abs[t/c-b ratio46[hi]]];rel=br/(t/c);
 <|"d"->d,"iterations"->its,"residual"->rr,"bracket_residual"->br,"relative_bracket_residual"->rel,"tolerance"->tol,"valid"->TrueQ[rel<=tol]|>];

(* g_{b,d}(0), the density of the centered b-sum.  The full centered
   characteristic function is used, not a Gaussian or binomial mixture. *)
gCenter46[b_Integer,d_?NumericQ,offset_?NumericQ,wp_Integer,step0_,cut0_]:=Module[{step=N[step0,wp],cut=N[cut0,wp],ww,logpsi,rr,g},
 ww=N[Range[-Round[cut/step],Round[cut/step]] step,wp];
 logpsi=b(-Log[1-I ww/Sqrt[b]]+Log[(1+Exp[-d+I d ww/Sqrt[b]])/(1+Exp[-d])])-I m46[d] ww Sqrt[b]-I ww offset/Sqrt[b];
 rr=Total[Re[Exp[logpsi]]];g=step rr/(2Pi Sqrt[b]);
 <|"g"->g,"log_g"->Log[g],"working_precision"->wp,"step"->step0,"cutoff"->cut0,"node_count"->Length[ww]|>];

(* Independent exact finite Gamma--binomial mixture in log-sum-exp form. *)
logFinite46[b_Integer,c_?NumericQ,t_?NumericQ,wp_Integer]:=Module[{cc=N[c,wp],tt=N[t,wp],x,ks,terms,mx,lp},
 x=tt/cc;ks=Range[0,b];terms=Table[If[x-k>0,LogGamma[b+1]-LogGamma[k+1]-LogGamma[b-k+1]+(b-1)Log[x-k],-Infinity],{k,ks}];
 mx=Max[terms];lp=mx+Log[Total[Exp[Select[terms,#>-Infinity&]-mx]]];
 -tt+(b-1)Log[cc]+lp-LogGamma[b]-b Log[1+Exp[-cc]]];
row46[b_Integer,c0_?NumericQ,cname_String,z_Integer,wp_Integer,step_,cut_,directQ_:False]:=Module[{c=N[c0,wp],t,sol,d,rho,eta,offset,gg,lf,ld,aa,logTail,logM,logPhysicalTail},
 t=N[b m46[c]+z Sqrt[b nu46[c]],wp];If[t<=0,Return[<|"b"->b,"c_label"->cname,"z"->z,"status"->"skipped: t<=0"|>]];
 sol=solveD46[b,c,t,wp];If[!TrueQ[sol["valid"]],Return[<|"b"->b,"c_label"->cname,"z"->z,"root_relative_bracket_residual"->ToString[N[sol["relative_bracket_residual"],20],InputForm],"status"->"root_failed"|>]];
 d=sol["d"];rho=d/c;eta=1-rho;offset=rho t-b m46[d];gg=gCenter46[b,d,offset,wp,step,cut];
 lf=b(Log[c/d]+Log[1+Exp[-d]]-Log[1+Exp[-c]])-eta t+Log[rho]+gg["log_g"];
 aa=N[cut/Sqrt[b],wp];logTail=-Log[Pi]-((b-2)/2)Log[1+aa^2]+Log[ArcTan[1/aa]];
 logM=-Log[rho]+Log[1+Exp[-d]]-Log[1+Exp[-c]];logPhysicalTail=Log[rho]+b logM-eta t+logTail;
 ld=If[directQ,logFinite46[b,c,t,wp],Missing["not_requested"]];
 <|"b"->b,"c_label"->cname,"c"->N[c,25],"z"->z,"t"->N[t,25],"d"->N[d,25],"rho"->N[rho,25],"eta"->N[eta,25],"center_offset_rho_t_minus_bm"->ToString[N[offset,20],InputForm],
 "logf_tilt_scale"->N[lf,25],"log_g_center"->N[gg["log_g"],25],"root_residual_string"->ToString[N[sol["residual"],20],InputForm],"root_bracket_residual_string"->ToString[N[sol["bracket_residual"],20],InputForm],"root_relative_bracket_residual_string"->ToString[N[sol["relative_bracket_residual"],20],InputForm],"root_relative_tolerance_string"->ToString[N[sol["tolerance"],20],InputForm],"root_valid"->sol["valid"],"root_iterations"->sol["iterations"],
 "working_precision"->wp,"fourier_step"->step,"fourier_cutoff"->cut,"fourier_nodes"->gg["node_count"],"CF_tail_a=cut/sqrt(b)"->N[aa,20],
 "log_center_density_tail_upper_bound"->N[logTail,20],"log_tail_over_gcenter"->N[logTail-gg["log_g"],20],"log_physical_density_tail_upper_bound"->N[logPhysicalTail,20],
 "independent_logsumexp_logf"->If[directQ,N[ld,25],Null],"independent_absolute_difference"->If[directQ,If[TrueQ[lf==ld],0,N[Abs[lf-ld],12]],Null],"status"->"completed"|>];

cs46={{N[4/5,constWP46],"0.8"},{mu46,"mu"},{N[9/5,constWP46],"1.8"}};zs46={-4,0,4};bs46={32,128,1024,65536};
(* The library guard exposes all kernel functions without launching scans. *)
If[!TrueQ[$blockKernel46LibraryMode],
base46=Flatten[Table[Module[{wp=If[b<=128,90,If[b<=1024,100,120]],cut=If[b==32,64,If[b==128,32,24]],step=1/8,cc=cs46[[j,1]],lab=cs46[[j,2]]},row46[b,cc,lab,z,wp,step,cut,b<=128]],{b,bs46},{j,Length[cs46]},{z,zs46}],2];
ref46=Table[Module[{a,brow},a=row46[65536,mu46,"mu",z,120,1/8,24,False];brow=row46[65536,mu46,"mu",z,240,1/16,32,False];
 <|"z"->z,"baseline"->a,"doubled"->brow,"logf_absolute_difference"->ToString[Abs[a["logf_tilt_scale"]-brow["logf_tilt_scale"]],InputForm],"logg_absolute_difference"->ToString[Abs[a["log_g_center"]-brow["log_g_center"]],InputForm]|>],{z,zs46}];
small46=Select[base46,#"b"<=128&];res46=<|"status"->"r046 stable tilted-scale block log-density evaluations; Fourier/trapezoid floating diagnostics, not an interval certificate",
 "wolfram"->$Version,"mu"->N[mu46,30],"formula"->"t=b m(c)+z sqrt(b nu(c)); solve t/c=b m(d)/d; log f=b[log(c/d)+log(1+e^-d)-log(1+e^-c)]-(1-d/c)t+log(d/c)+log g_(b,d)(rho t-b m(d)); exact roots make the offset zero",
 "centered_CF"->"psi_d(w)=exp(-I w m(d))/(1-I w)*(1+exp(-d+I d w))/(1+exp(-d)); integrate psi_d(w/sqrt(b))^b exp(-I w offset/sqrt(b))/(2 pi sqrt(b))",
 "CF_tail_bound"->"beyond |w|>cut in normalized Fourier variable: tail density <= Pi^-1 (1+(cut/sqrt(b))^2)^(-(b-2)/2) ArcTan[sqrt(b)/cut). Physical tilted-density tail multiplies by rho Exp[b log M_c(eta)-eta t]. This does not certify the trapezoid discretization error.",
 "base_grid"->base46,"small_b_logsumexp_check"-><|"rows"->Length[small46],"max_absolute_logf_difference"->N[Max[small46[[All,"independent_absolute_difference"]]],12],"cutoffs"-><|"b32"->64,"b128"->32|>|>,
 "b65536_precision_cutoff_step_doubling"->ref46,"run_status"->"completed"|>;
Export[out46,res46,"RawJSON"];Print[<|"rows"->Length[base46],"small_max_diff"->res46["small_b_logsumexp_check","max_absolute_logf_difference"],"large_checks"->ref46[[All,{"z","logf_absolute_difference"}]]|>];
];
