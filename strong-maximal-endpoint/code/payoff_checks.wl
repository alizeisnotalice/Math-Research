Clear[x,c,theta,r];
M[c_,theta_]:=(1+Exp[-(1-theta)c])/((1-theta)(1+Exp[-c]));
stationaryResidual=FullSimplify[D[Log[M[c,theta]],c]-(1/(1+Exp[c])-(1-theta)/(1+Exp[(1-theta)c])),Assumptions->{c>0,0<theta<1}];
hazardResiduals=Table[With[{pp=Sum[x^j/Factorial[j],{j,0,b-1}]},Expand[(b-1)pp-x D[pp,x]-Sum[(b-1-j)x^j/Factorial[j],{j,0,b-1}]]],{b,1,8}];
chernoffHalf=FullSimplify[4Exp[-3]<1/2];mu=1+ProductLog[Exp[-1]];
penaltyRoots=Table[With[{rr=nn+ProductLog[-nn Exp[-nn]]},<|"N"->nn,"positive_penalty_root"->ToString[N[rr,50],InputForm],"N_minus_root"->ToString[N[-ProductLog[-nn Exp[-nn]],50],InputForm]|>],{nn,{2,16,64,128}}];
frontRoots=Table[With[{a=1+2/nn},With[{qq=-a-ProductLog[-1,-Exp[-a]]},<|"N"->nn,"first_stage_front"->ToString[N[qq,50],InputForm],"equation_residual"->With[{qn=N[qq,110]},ToString[N[Exp[qn]-qn-1-2/nn,50],InputForm]]|>]],{nn,{16,32,64,128}}];
mgfRows=Table[Module[{cc},cc=cv/.FindRoot[Exp[(1-th)cv]-(1-th)Exp[cv]+th==0,{cv,N[mu,70]},WorkingPrecision->60,AccuracyGoal->40,PrecisionGoal->40];<|"theta"->ToString[th,InputForm],"global_max_c"->ToString[N[cc,40],InputForm],"Mstar"->ToString[N[M[cc,th],40],InputForm]|>],{th,{1/20,1/10,1/4,1/2}}];
result=<|"wolfram"->$Version,"stationarity_residual"->ToString[stationaryResidual,InputForm],"Gamma_hazard_polynomial_residuals"->hazardResiduals,"four_e_minus3_less_half"->chernoffHalf,"penalty_roots"->penaltyRoots,"first_front_roots"->frontRoots,"positive_mgf_maxima"->mgfRows|>;


Clear[aa,vf,rr,cc];dd=aa[rr]+Exp[-cc]aa[rr-cc]-(1+Exp[-cc])vf[rr-cc];dResidual=FullSimplify[(D[dd,cc]/.Derivative[1][aa][zz_]:>vf[zz]-aa[zz])-(1+Exp[-cc])Derivative[1][vf][rr-cc]];result=Join[result,<|"D_derivative_residual"->ToString[dResidual,InputForm],"half_tilt_Mstar_residual"->FullSimplify[M[2Log[1+Sqrt[2]],1/2]-(1+Sqrt[2])]|>];structureResult=result;
mu=1+ProductLog[Exp[-1]];c0=N[mu,70];p0=1/(1+Exp[c0]);
parts[a_Integer,b_Integer,r_?NumericQ,nn_Integer]:=Module[{ks=Range[0,b],bw,x,u,q,call,rr},If[r<=0,Return[{0,1,a+b c0 p0-r}]];bw=N[Binomial[b,#]&/@ks,70]Exp[-c0 ks]/(1+Exp[-c0])^b;x=r-c0 ks;u=Clip[x,{0,Infinity}];q=Total[bw GammaRegularized[a,u]];call=Total[bw(a GammaRegularized[a+1,u]-x GammaRegularized[a,u])];rr=Exp[-r]Total[N[Binomial[b,#]&/@ks,70]u^a]/((1+Exp[-c0])^b Gamma[a+1]);{rr,q,call}];
v[a_Integer,b_Integer,r_?NumericQ,nn_Integer]:=If[a==0&&b==0,If[r>0,Exp[-r],1-r/nn],With[{z=parts[a,b,r,nn]},z[[1]]+z[[2]]+z[[3]]/nn]];
vp[a_Integer,b_Integer,r_?NumericQ,nn_Integer]:=If[a==0&&b==0,If[r>0,-Exp[-r],-1/nn],With[{z=parts[a,b,r,nn]},-z[[1]]-z[[2]]/nn]];
density[t_Integer,r_?NumericQ]:=If[r<=0,0,Exp[-r]Total[N[Binomial[t,#]&/@Range[0,t],70]Clip[r-c0 Range[0,t],{0,Infinity}]^(t-1)]/((1+Exp[-c0])^t Gamma[t])];
fprime[t_Integer,r_?NumericQ]:=Module[{x=Clip[r-c0 Range[0,t],{0,Infinity}],bw=N[Binomial[t,#]&/@Range[0,t],70]},density[t,r](-1+(t-1)Total[bw x^(t-2)]/Total[bw x^(t-1)])];
case[nn_,tt_,rr_]:=Module[{aa,cc=c0,lc=0,hc=rr+c0,d,diff,sl,proposal,base,b1,b2,b3,exactgain,taylor,src,pr,ff,fp,z,iter=0,kk},aa=v[tt,tt-1,rr,nn];base=v[tt,tt,rr,nn];Do[iter=k;diff=aa+Exp[-cc]v[tt,tt-1,rr-cc,nn]-(1+Exp[-cc])v[tt-1,tt-1,rr-cc,nn];If[diff>0,lc=cc,hc=cc];sl=(1+Exp[-cc])vp[tt-1,tt-1,rr-cc,nn];If[Abs[diff]<10^-35 Max[Abs[sl],10^-30],Break[]];proposal=cc-diff/sl;If[proposal<=lc||proposal>=hc,proposal=(lc+hc)/2];If[Abs[proposal-cc]<10^-40,Break[]];cc=proposal,{k,1,30}];
kk=(aa+Exp[-cc]v[tt,tt-1,rr-cc,nn])/(1+Exp[-cc]);exactgain=kk-base;
b1=p0(base-v[tt-1,tt-1,rr-c0,nn]);b2=p0(2p0-1)(base-v[tt-1,tt-1,rr-c0,nn])+p0 vp[tt-1,tt-1,rr-c0,nn];taylor=-b1^2/(2b2);
pr=parts[tt,tt,rr,nn];ff=density[tt,rr];fp=fprime[tt,rr];b3=ff-pr[[1]]-(1-1/nn)fp;src=mu^2 p0/8 (pr[[1]]-(1-1/nn)ff)^2/(pr[[1]]+pr[[2]]/nn);
<|"N"->nn,"future_t"->tt,"r"->N[rr,20],"R_over_f"->N[pr[[1]]/ff,20],"v2_over_f"->N[(pr[[1]]-(1-1/nn)ff)/ff,20],"v3_over_f"->N[b3/ff,20],"best_c"->N[cc,20],"root_iterations"->iter,"root_residual"->ToString[N[aa+Exp[-cc]v[tt,tt-1,rr-cc,nn]-(1+Exp[-cc])v[tt-1,tt-1,rr-cc,nn],30],InputForm],"actual_max_gain"->N[exactgain,20],"exact_parameter_Taylor_gain"->N[taylor,20],"actual_over_Taylor"->N[exactgain/taylor,20],"spatial_source_candidate"->N[src,20]|>];
crit=Table[With[{z=zv/.FindRoot[Erfc[zv/Sqrt[2]]/2-Sqrt[mu]zv/Sqrt[nn]==0,{zv,1},WorkingPrecision->70]},case[nn,nn,nn mu-Sqrt[nn mu]z]],{nn,{16,64,256}}];
moderate=Table[case[256,32,32mu+zz Sqrt[32mu]],{zz,{-4,0,4}}];
result=<|"wolfram"->$Version,"status"->"70-digit explicit Gamma-binomial expectation, one-step global root; no interval certificate","critical_rows"->crit,"moderate_rows"->moderate|>;


localResult=result;result=<|"structure"->structureResult,"local_gain"->localResult|>;base=DirectoryName[DirectoryName[$InputFileName]];Export[FileNameJoin[{base,"rounds","r033","symbolic_checks.json"}],result,"RawJSON"];Print["r033 structure and one-step checks completed"];
