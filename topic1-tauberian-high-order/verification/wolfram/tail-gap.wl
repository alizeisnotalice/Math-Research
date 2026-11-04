Print["VERSION ",$Version];
(* h=4^n >=16, D,w>=0, B=8(hD+D+2)+w; 0<v, Bv<1/16, 4h(D+1)v<1, 1/2<=r<=1-(B-2D)v. *)
B=8(h D+D+2)+w; gap=(1-r)/2-v(1-B v)/(1+v);
certificate=(1-r-(B-2D)v)/2+v^2(1+B)/(1+v)+v(15h D/4+3D+7+w/2);
Print["gap_certificate_identity=",InputForm[FullSimplify[Together[gap-h D v/4-certificate]==0]]];
Print["cap_slack_positive=",InputForm[Resolve[ForAll[{h,D,v,r},Implies[h>=16&&D>=0&&v>0&&4h(D+1)v<1&&r>=1/2,r/4>h D v/4]],Reals]]];
Print["radius_series=",InputForm[Normal[Series[(1-B v-v^3)/((1+v)(1-D v)^(1/n)),{v,0,1}]]]];
