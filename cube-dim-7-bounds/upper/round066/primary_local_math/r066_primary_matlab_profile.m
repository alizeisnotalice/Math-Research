% New first-profile antiderivative and target-Fubini primitive identities.
% Preserve the existing shared workspace and path; use a fresh variable prefix.
r066_started=tic;
syms r066_R positive
syms r066_u real
syms r066_a positive
r066_Q=-log(1-r066_u/(2*r066_R));
r066_profile_primitive=r066_R*r066_Q^2;
r066_identity_profile=simplify(diff(r066_profile_primitive,r066_u)-r066_Q/(1-r066_u/(2*r066_R)));
r066_tail_primitive=1-r066_a+r066_a*log(r066_a);
r066_identity_mean=simplify(diff(r066_tail_primitive,r066_a)-log(r066_a));
r066_mu=3*subs(r066_tail_primitive,r066_a,sym(5)/6);
r066_identity_mu=simplify(r066_mu-(sym(1)/2-sym(5)/2*log(sym(6)/5)));
assert(isAlways(r066_identity_profile==0));
assert(isAlways(r066_identity_mean==0));
assert(isAlways(r066_identity_mu==0));
r066_receipt=struct('scope','New symbolic identities for Jphi tail primitive and target mean; no full operator proof inferred',...
 'profile_antiderivative_identity',char(r066_identity_profile),...
 'mean_tail_identity',char(r066_identity_mean),'relative_mean_identity',char(r066_identity_mu),...
 'mu_exact',char(r066_mu),'mu_50digits',char(vpa(r066_mu,50)),...
 'matlab_version',version,'wall_seconds',toc(r066_started),...
 'time_utc',char(datetime('now','TimeZone','UTC')),'workspace_cleared',false,'path_reset',false,...
 'domain_note','Real positivity and capture/Fubini hypotheses come from the analytic proof, not a symbolic simplify call.');
r066_fid=fopen('/Users/zhengzhihao/Desktop/T/cube-v3-work/r066_primary_matlab_receipt.json','w');assert(r066_fid~=-1);fprintf(r066_fid,'%s\n',jsonencode(r066_receipt));fclose(r066_fid);
disp(r066_receipt)
