base=DirectoryName[$InputFileName];
run[task_]:=Module[{nn=task[[1]],wp=task[[2]],mu,p,m,b,t,ww,lg,gam,den,num,sc,su,cp,out},
 mu=N[1+ProductLog[Exp[-1]],wp+20];p=1/(1+Exp[mu]);m=mu nn;b=N[nn^(1/4),wp+20];t=N[nn^(-3/4),wp+20];
 ww=Table[N[Binomial[nn,j]p^j(1-p)^(nn-j),wp+10],{j,0,nn}];lg=N[LogGamma[nn],wp+20];
 gam[x_]:=If[x<=0,0,Exp[(nn-1)Log[x]-x-lg]];
 den=Sum[ww[[j+1]]gam[m-mu j],{j,0,nn}];num=Sum[ww[[j+1]]Exp[-(nn-j)t]gam[m-mu j-(nn-j)t],{j,0,nn}];
 sc=num/den;su=(p+(1-p)Exp[-t])^nn;
 cp=Sum[If[m+2b-mu j<=0,0,ww[[j+1]]Gamma[nn,Max[0,m-2b-mu j],m+2b-mu j]/Gamma[nn]],{j,0,nn}];
 <|"n"->nn,"precision_digits"->wp,"conditional_survival"->ToString[N[sc,wp],InputForm],"unconditional_survival"->ToString[N[su,wp],InputForm],"central_probability"->ToString[N[cp,wp],InputForm],"tv_lower_expression"->ToString[N[1-sc-su-cp,wp],InputForm],"binomial_normalization_residual"->ToString[N[Total[ww]-1,wp],InputForm]|>
];
LaunchKernels[8];DistributeDefinitions[run];res=ParallelMap[run,Flatten[Table[{nn,wp},{nn,{16,64,256,1024}},{wp,{50,80}}],1]];
Export[FileNameJoin[{base,"rank_tv_wolfram.json"}],<|"status"->"complete","kernels"->Length[Kernels[]],"scope"->"High precision evaluation of finite positive formulas; not interval certification.","results"->res|>,"RawJSON"];
Print[<|"status"->"complete","cases"->Length[res]|>];CloseKernels[];Quit[];
