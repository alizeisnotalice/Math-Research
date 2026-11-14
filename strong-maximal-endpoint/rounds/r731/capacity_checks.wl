primal=Maximize[{x^2 y,x+y==2&&x>=0&&y>=0},{x,y}];
Print["star_geometric_primal=",primal];
value=First[primal]^(1/3);
kl=(2/3)Log[4/3]+(1/3)Log[2/3];
Print["primal_log_dual_residual=",FullSimplify[Log[value]-kl]];
ff=(2/3)Log[4/3]+(t/3)Log[4t/3]+((1-t)/3)Log[4(1-t)/3]-t Log[2t]-(1-t)Log[2(1-t)];
Print["hessian_identity_residual=",FullSimplify[D[ff,{t,2}]+2/(3t(1-t)),Assumptions->0<t<1]];
Print["symmetric_gradient=",FullSimplify[D[ff,t]/.t->1/2]];
Quit[];
