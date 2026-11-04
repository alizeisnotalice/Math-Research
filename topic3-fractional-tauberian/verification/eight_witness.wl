ClearAll[a,t,w,x,y];
$Assumptions=Element[{a,t,w,x,y},Reals]&&1/Sqrt[2]<a<1&&a<w<=1;
c=2a^2-1;A0=(1/a-a)/2;B0=(a-Sqrt[c])/2;aa=(w+1/w)/2-a;bb=(1+w^2-2a w)/(4a);
ws=2a^2/(a+Sqrt[c]);sigma=(ws-a)(1-a ws)/(ws(1-a^2)(1-a));
res=FullSimplify/@{bb/aa-w/(2a),ws-2a B0/A0,(a-2B0)^2-c,
 (bb-(1/w-w)/4)-(w-a)(1+w^2)/(4a w)};
If[res =!= {0,0,0,0},Print[res];Exit[1]];
v2=x y;v3=x^2+y^2-x y-4t^2;v4=(x-2t)(y-2t);
dual=Expand/@{3+v4-2(1-t)(x+y),1+v2+v3+v4-2(1-t)(x+y)};
scoef=FullSimplify/@Table[SeriesCoefficient[sigma/.a->1-t,{t,0,j}],{j,0,12}];
wcoef=FullSimplify/@Table[SeriesCoefficient[ws/.a->1-t,{t,0,j}],{j,0,14}];
polynomial=FullSimplify[(1-a^2)ws^2-4a^3 ws+4a^4];If[polynomial=!=0,Exit[2]];
root=Root[4#^3-#-1&,1];
out=<|"version"->$Version,"identities"->res,"w_implicit_residual"->polynomial,
 "alpha_transition"->ToString[N[root,70],InputForm],
 "sigma_coefficients"->(ToString[#,InputForm]&/@scoef),"w_coefficients"->(ToString[#,InputForm]&/@wcoef),
 "dual_polynomials"->(ToString[#,InputForm]&/@dual),
 "scope"->"exact symbolic checks for finite-witness rectangle model and four-axis capacity"|>;
Export[FileNameJoin[{DirectoryName[DirectoryName[$InputFileName]],"results","eight_witness_wolfram.json"}],out,"JSON"];Print[out];

ClearAll[sg];at=1-t;qt=(4-Sqrt[7])/6;pt=2qt-qt^2;gam=FullSimplify[2qt-2qt^2];rad0=(1-at^2)/(at+Sqrt[at^2-3(1-at^2)/4]);rr=rad0(1-sg t);ee=Sqrt[1+(1/4+4pt)rr^2]-rr;ee2=ee+2(1-qt)rr;dlow=ee(ee-at)-pt rr^2;
ld=<|"gamma"->ToString[gam,InputForm],"mass_residual"->FullSimplify[(ee+3rr/2)(ee+rr/2)-4pt rr^2-1],"low_slack_leading"->ToString[FullSimplify[Limit[dlow/t^2,t->0,Direction->"FromAbove"]],InputForm],"radial_excess_leading"->ToString[FullSimplify[Limit[(rr-(1/at-at)/2)/t^2,t->0,Direction->"FromAbove"]],InputForm]|>;
ld["finite_signs"]=Table[With[{dv=(dlow/.{sg->7/20,t->u}),rv=((rr-(1/at-at)/2)/.{sg->7/20,t->u})},{ToString[u,InputForm],Sign[RootReduce[dv]],Sign[RootReduce[rv]]}],{u,{1/1000,1/10000,1/100000}}];
If[ld["mass_residual"]=!=0 || ld["finite_signs"][[2,2]]=!=1,Exit[3]];out["L_corner"]=ld;out["cauchy_bound_M"]=4;out["cauchy_radius"]="1/8";
Export[FileNameJoin[{DirectoryName[DirectoryName[$InputFileName]],"results","eight_witness_wolfram.json"}],out,"JSON"];Print[ld];
