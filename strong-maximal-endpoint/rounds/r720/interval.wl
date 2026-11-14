folder=DirectoryName[$InputFileName];u={1,3,4,3,1};
wts={{1/4,5/4,1,5/4,1/4},{1/4,19/20,7/10,17/10,1},{1,7/5,2/5,7/5,1}};
exactLogExpand[e_]:=Expand[e/.{Log[r_Rational]:>Total[(Last[#]Log[First[#]]&)/@FactorInteger[r]],Log[r_Integer]:>Total[(Last[#]Log[First[#]]&)/@FactorInteger[r]]}];
res=Table[Module[{w=wts[[i]],ex,z,d},ex=20u w;z=Times@@MapThread[Power,{w,ex}];d=Total[u w Log[w]]/12;<|"case"->i,"integer_exponents"->And@@(IntegerQ/@ex),"exp240KL"->ToString[z,InputForm],"KL100"->ToString[N[d,100],InputForm],"exact_log_identity"->(exactLogExpand[Log[z]-240d]==0)|>],{i,3}];
z=ToExpression[# ["exp240KL"]]&/@res;ratio={z[[2]]/z[[1]],z[[3]]/z[[2]]};
ell0=(6Log[7/5]+4Log[2/5])/10;ellArm=3Log[7/5]/4;pin=ellArm-ell0-2(3/12)(7/10)^2;
out=<|"version"->$Version,"results"->res,"successive_ratio_floors"->(Floor/@ratio),"strict_order"->(ratio[[1]]>1&&ratio[[2]]>1),"responsible_LF_test_volume"->ToString[24/5,InputForm],"max_KL_test_volume"->ToString[30/7,InputForm],"price_loss_ratio"->ToString[(24/5)/(30/7),InputForm],"weight_TV"->ToString[Total[Abs[({10,1,1}-{4,4,4})/12]]/2,InputForm],"whole_overlap_Pinsker_margin100"->ToString[N[pin,100],InputForm]|>;
Export[FileNameJoin[{folder,"wolfram.json"}],out,"RawJSON"];Print[<|"integer_log_identities"->And@@(#["exact_log_identity"]&/@res),"strict_order"->out["strict_order"],"ratio_floors"->out["successive_ratio_floors"]|>];Quit[];
