function result=sm_laplace_checks
% Exact Laplace minimizer; quantile optimization remains a numerical diagnostic.
root=fileparts(fileparts(mfilename('fullpath')));
B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B); sigma=sqrt(mu);
laplace=cell(0,1);
for lambda=[.5 .2 .1 .05 .02 .01 .005 .002 .001]
 station=@(c)exp(c).*(expm1(lambda*c)/lambda-1)-1;
 c=fzero(station,[0 mu],optimset('TolX',1e-13));
 derivative=@(s)lambda*(exp(s).*(expm1(lambda*s)/lambda-1)-1) ...
     ./((1+exp(s)).*(1+exp((1+lambda)*s)));
 gain=integral(derivative,mu,c,'AbsTol',1e-18,'RelTol',1e-10);
 gainDirect=logM(c,lambda)-logM(mu,lambda);
 row=struct('lambda',lambda,'optimal_c',c,'stationarity_residual',station(c), ...
  'linear_shift_ratio',(c-mu)/lambda,'cubic_gain_ratio',gain/lambda^3, ...
  'integral_direct_difference',gain-gainDirect);
 laplace{end+1}=row; %#ok<AGROW>
 fprintf('lambda=%g c=%.12g shift/lambda=%.9g gain/lambda^3=%.9g\n', ...
  lambda,c,row.linear_shift_ratio,row.cubic_gain_ratio);
end
quantiles=cell(0,1);
for n=[16 32 64 128 256 512 1024 2048 4096]
 beta=sigma*sqrt(log(n)/n); q0=quantile(n,mu,beta,mu);
 [c,negq,exitflag]=fminbnd(@(c)-quantile(n,c,beta,mu),.3,2.3,optimset('TolX',1e-8));
 [~,~,residual]=sm_tail(n,exp(-c),-negq);
 quantiles{end+1}=struct('n',n,'beta',beta,'candidate_c',c, ...
  'fixed_mu_quantile',q0,'candidate_quantile',-negq,'gain',-negq-q0, ...
  'gain_over_log_n',(-negq-q0)/log(n),'candidate_shift_over_rho', ...
  (c-mu)/sqrt(log(n)/n),'quantile_residual',residual-beta,'exitflag',exitflag); %#ok<AGROW>
 fprintf('n=%4d gain=%.10g gain/logn=%.9g flag=%d\n',n,-negq-q0,(-negq-q0)/log(n),exitflag);
end
result=struct('matlab_version',version,'precision','double; not directed certification', ...
 'parameters',struct('quantile_c_bounds',[.3 2.3],'optimizer_tol_x',1e-8, ...
 'quantile_bisections',55,'laplace_root_tol_x',1e-13,'randomness','none'), ...
 'predicted_laplace_shift_limit',-mu/2,'proved_cubic_gain_limit',-mu*(mu-1)/8, ...
 'unproved_quantile_gain_coefficient',(mu-1)/8,'laplace',{laplace},'quantiles',{quantiles});
fid=fopen(fullfile(root,'rounds','r015','checks.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
function value=logM(c,lambda)
value=log1p(exp(-(1+lambda)*c))-log1p(exp(-c))-log1p(lambda);
end
function q=quantile(n,c,beta,mu)
lo=0; hi=n*mu+12*sqrt(n);
for j=1:55
 mid=(lo+hi)/2; [~,~,failure]=sm_tail(n,exp(-c),mid);
 if failure<beta,lo=mid;else,hi=mid;end
end
q=(lo+hi)/2;
end
