function result=sm_edgeworth_ramp_checks
% Reuse archived soft-root computations; no expensive mixture rerun.
root=fileparts(fileparts(mfilename('fullpath')));
old=jsondecode(fileread(fullfile(root,'rounds','r017','checks.json')));
B=fzero(@(b)log(b)-1-1/b,[3 4]);mu=log(B);sigma=sqrt(mu);
a=(3-mu)/6;b=(5*mu-6)/6;rows=cell(0,1);
for r=old.rows'
 n=r.n;L=log(n);fun=@(z)log(erfcx(z/sqrt(2))/2)-z.^2/2-log(sigma)-log(z)+L/2;
 lo=realmin;hi=sqrt(L)+2;
 for j=1:70
  mid=(lo+hi)/2;if fun(mid)>0,lo=mid;else,hi=mid;end
 end
 z=(lo+hi)/2;dn=sigma*sqrt(n)*z;
 correction=r.fixed_mu_q-n*mu+dn;
 asymptoticLog=n*mu-dn+a*z^2+b-.5*log(pi*n/2);
 score=r.fixed_mu_q-r.log_I_Stirling;
 row=struct('n',n,'z',z,'fixed_soft_constant_residual',correction-a*z^2, ...
  'constant_target',b,'fixed_score_over_family_asymptotic',exp(score-asymptoticLog), ...
  'fixed_score_log_n_coefficient',(score-n*mu+dn)/L, ...
  'shifted_score_log_n_coefficient',(r.shifted_q-r.log_I_Stirling-n*mu+dn)/L);
 rows{end+1}=row; %#ok<AGROW>
 fprintf('n=%g constant=%.10g target=%.10g family_ratio=%.10g\n', ...
 n,row.fixed_soft_constant_residual,b,row.fixed_score_over_family_asymptotic);
end
result=struct('matlab_version',version,'status','double-precision diagnostics derived from archived r017 data', ...
 'fixed_log_n_target',-mu/6,'shifted_log_n_target',-(mu+3)/24, ...
 'family_prefactor',sqrt(2/pi)*exp(b),'rows',{rows});
fid=fopen(fullfile(root,'rounds','r018','checks.json'),'w');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
