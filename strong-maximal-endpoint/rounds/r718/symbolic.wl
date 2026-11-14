folder=DirectoryName[$InputFileName];
y0=(aa+1)ell/(cc rt);assume=aa>=1&&ell>0&&cc>0&&rt>0&&eta>y0;
integral=Integrate[Exp[aa ell-cc rt y],{y,y0,eta},Assumptions->assume];
identity=FullSimplify[integral-(Exp[-ell]-Exp[aa ell-cc rt eta])/(cc rt),Assumptions->assume];
W=1/(rt ell);h=W/nn;rho=h rt/(2cc ell);
checks=<|"tail_integral_identity"->identity,"grid_threshold_identity"->FullSimplify[rho-1/(2cc ell^2 nn)],"large_branch_log_power_lower"->Limit[Log[2cc ell^4]/Log[ell],ell->Infinity,Assumptions->cc>0],"large_branch_log_power_upper"->Limit[Log[2cc ell^2(ell^2+1)]/Log[ell],ell->Infinity,Assumptions->cc>0],"Markov_retained_fraction"->FullSimplify[1-(cc rr ell/rt)/(2cc rr ell/rt),Assumptions->rr>0&&ell>0&&rt>0&&cc>0],"direct_score_ratio"->FullSimplify[(uu/2)/ss/(uu/mm)-1/(2(ss/mm))],"near_median_scaled_width_limit"->Limit[1/ell,ell->Infinity]|>;
checks["Markov_retained_fraction"]=ToString[checks["Markov_retained_fraction"],InputForm];
Export[FileNameJoin[{folder,"wolfram.json"}],<|"version"->$Version,"checks"->checks,"scope"->"Symbolic identities only; no certification of unknown Gaussian constants or starting dimension"|>,"RawJSON"];Print[checks];Quit[];
