Print["VERSION ",$Version," PID ",$ProcessID];
(* R>=1, candidate radius r>0, center distance c>=0, and a unit-radius root at origin. *)
Print["radius_violation=",InputForm[Reduce[R>=1 && r>0 && c>=0 && c<=R+r/2 && c>Max[1,r] && r>=2R,{R,r,c},Reals]]];
Print["center_violation=",InputForm[Reduce[R>=1 && r>0 && c>=0 && c<=R+r/2 && c>Max[1,r] && c>=2R,{R,r,c},Reals]]];
Print["coverage_separation_violation=",InputForm[Reduce[r>0 && dist>=0 && dist<=1/2 && dist>Max[1,r],{r,dist},Reals]]];
Print["radius_m_violation=",InputForm[Reduce[R>=1 && m>=8R+8 && 1/2+m/4<=2R+1,{R,m},Reals]]];
Print["ceiling_radius_lower_slack=",InputForm[FullSimplify[(1/2+(8R+8)/4)-(2R+1),R>=1]]];
(* Strictness matters; if the center condition is weakened to >=, r=c=2R is feasible. *)
Print["nonstrict_boundary_witness=",InputForm[FullSimplify[{c<=R+r/2,c>=Max[1,r],r>=2R}/.{R->1,r->2,c->2}]]];
