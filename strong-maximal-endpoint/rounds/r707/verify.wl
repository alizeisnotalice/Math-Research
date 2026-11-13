mu=1+ProductLog[1/E];a0=(3-mu)/6;d0=(mu-1-Log[mu])/2;aq=(mu+3)/8-Log[2 Pi mu];
rOverZ=x-x^2+3 x^3-15 x^4+105 x^5;hh=1/(1+rOverZ);
soft=Series[hh*(a*(1/x-1)+u*(1+(1-1/x)*rOverZ/2)),{x,0,2}];
ratio=Series[hh^2/x,{x,0,2}];
combined=FullSimplify[(7+5u)/8-u-Log[2 Pi Sqrt[u]]+(u-1-Log[u])/2-((u+3)/8-Log[2 Pi u]),Assumptions->u>0];
large=Table[Module[{n=10^j,zz,z,rr,h,nn,tt,gain,slab,reduced,expanded},
 z=zz/.FindRoot[Log[Erfc[zz/Sqrt[2]]/2]==Log[Sqrt[mu/n]zz],{zz,Sqrt[Log[n]]},WorkingPrecision->120];
 rr=(Erfc[z/Sqrt[2]]/2)/(Exp[-z^2/2]/Sqrt[2 Pi]);h=1/(1+rr/z);nn=z/rr;
 tt=h*(a0*(z^2-1)+mu*(1+(1-z^2)*rr/(2z)));gain=(mu-1)*z^2*h^2/8;
 slab=(mu*z^2/2-Log[n]/2-Log[nn]-Log[2 Pi]/2)/(1+nn)-1/nn;
 reduced=tt+gain+slab+d0/Log[n];expanded=(9-mu)*z^2/24+(13mu-15)/12+(-4Log[z]+aq)/z^2;
 <|"log10_n"->j,"z"->ToString[z,InputForm],"scaled_resummed_remainder"->ToString[N[z^2*(reduced-expanded),80],InputForm]|>
 ],{j,{6,10,20,50,100,200,500,1000}}];
cdf=Flatten[Table[Module[{xx=N[n+t Sqrt[n],90],val},
 val=TimeConstrained[GammaRegularized[N[n,80],0,xx],20,$Failed];
 <|"n"->n,"standardized_x"->t,"x"->ToString[xx,InputForm],"cdf"->ToString[val,InputForm]|>
 ],{n,{4194304,16777216,67108864}},{t,{-4,-7/2,-3,-2,0,4}}]];
targets=Import[FileNameJoin[{DirectoryName[$InputFileName],"gamma_certificates.json"}],"RawJSON"]["rows"];
targeted=Flatten[Table[Module[{xx=rec["x_numerator"]/rec["x_denominator"],n=rec["n"],value},value=TimeConstrained[GammaRegularized[N[n,prec],0,N[xx,prec+10]],20,$Failed];<|"n"->n,"shift"->rec["shift"],"precision"->prec,"cdf"->ToString[value,InputForm]|>],{rec,Select[targets,#["shift"]==0&]},{prec,{80,120}}]];
out=<|"version"->$Version,"targeted_cdf_checks"->targeted,"soft_series"->ToString[soft,InputForm],"soft_one_over_z2_identity"->ToString[FullSimplify[3a0+mu/2],InputForm],
 "ratio_series"->ToString[ratio,InputForm],"combined_identity"->ToString[combined,InputForm],
 "A_Q"->ToString[N[aq,90],InputForm],"large_dimension_reduced_checks"->large,"gamma_cdf_checks"->cdf|>;
Export[FileNameJoin[{DirectoryName[$InputFileName],"wolfram.json"}],out,"RawJSON"];Print[out];Quit[];
