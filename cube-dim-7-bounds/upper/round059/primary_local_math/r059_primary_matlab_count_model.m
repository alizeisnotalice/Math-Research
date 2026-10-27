% New R059 symbolic count calibration and ideal-model gradient diagnostics.
% Existing shared MATLAB workspace/path are preserved. No reset or clear.
r059_clock=tic;
syms r059_q r059_D r059_c r059_n r059_r real
r059_ell=(r059_q+r059_D+r059_c*r059_q)/(2*(r059_q+r059_D));
r059_identity_ell=simplify(diff(r059_ell,r059_q)-r059_c*r059_D/(2*(r059_q+r059_D)^2));
r059_count=r059_n+(r059_n-1)*r059_c*r059_q/r059_r;
r059_ratio=r059_count/(2*r059_n*r059_ell);
r059_expected=(r059_q+r059_D)*(r059_r+r059_c*r059_q)/(r059_r*(r059_q+r059_D+r059_c*r059_q))-(r059_q+r059_D)*r059_c*r059_q/(r059_n*r059_r*(r059_q+r059_D+r059_c*r059_q));
r059_identity_count=simplify(r059_ratio-r059_expected);
assert(isAlways(r059_identity_ell==0)&&isAlways(r059_identity_count==0));
syms r059_d r059_a r059_b real
r059_qmin=1-r059_b*r059_d;r059_qmax=1-r059_a*r059_d;
r059_Rmin=r059_qmin+r059_d;
r059_phi=-log(1-r059_d/(2*r059_Rmin));
r059_mu=r059_Rmin/r059_qmin*r059_phi^2;
r059_cap=r059_a*r059_d/(2*r059_qmax);
r059_model=log((r059_qmax+r059_d)/r059_Rmin)*r059_qmin/r059_qmax-2*r059_mu;
r059_tuple=[sym(391)/1905,sym(283)/2032,sym(63)/64];
r059_vars=[r059_d,r059_a,r059_b];
r059_output=struct('scope','Two exact symbolic count identities and ideal-model local gradient; no operator bound or optimum claim',...
 'count_identity',char(r059_identity_count),'ell_derivative_identity',char(r059_identity_ell),...
 'tuple',{{'391/1905','283/2032','63/64'}},...
 'model_value',char(vpa(subs(r059_model,r059_vars,r059_tuple),50)),...
 'mu',char(vpa(subs(r059_mu,r059_vars,r059_tuple),50)),...
 'cap_minus_mu',char(vpa(subs(r059_cap-r059_mu,r059_vars,r059_tuple),50)),...
 'model_gradient_d_a_b',{{char(vpa(subs(diff(r059_model,r059_d),r059_vars,r059_tuple),50)),char(vpa(subs(diff(r059_model,r059_a),r059_vars,r059_tuple),50)),char(vpa(subs(diff(r059_model,r059_b),r059_vars,r059_tuple),50))}},...
 'matlab_version',version,'wall_seconds',toc(r059_clock),'path_reset',false,'workspace_cleared',false);
r059_fid=fopen('/Users/zhengzhihao/Desktop/T/cube-v3-work/r059_primary_matlab_receipt.json','w');assert(r059_fid~=-1);fprintf(r059_fid,'%s\n',jsonencode(r059_output));fclose(r059_fid);
disp(r059_output)
