function result=sm_independent_log_checks
B=fzero(@(b)log(b)-1-1/b,[3 4]);mu=log(B);sig=sqrt(mu);costRows=cell(0,1);
for lambda=[.2 .1 .05 .02 .01 .005 .002 .001]
 objective=@(c)mu-meanAt(c)+lambda*(varianceAt(c)-mu)/2;
 [c,value,flag]=fminbnd(objective,0,20,optimset('TolX',1e-10));
 costRows{end+1}=struct('lambda',lambda,'candidate_c',c, ...
 'mean_variance_cost_over_lambda_squared',value/lambda^2, ...
 'shift_over_lambda',(c-mu)/lambda,'exitflag',flag); %#ok<AGROW>
 fprintf('lambda=%g minH/lambda²=%.10g shift/lambda=%.10g\n',lambda,value/lambda^2,(c-mu)/lambda);
end
rows=cell(0,1);
for n=[2 4 8 16 32 64 128 256]
 ell=log(n);rho=sqrt(ell/n);c1=mu-sig*rho/2;c2=mu+sig*rho/3;n1=floor(n/2);n2=n-n1;
 p1=1/(1+exp(c1));p2=1/(1+exp(c2));k=(0:n1)';j=0:n2;
 w1=binopdf(k,n1,p1);w2=binopdf(j,n2,p2);logW=log(w1)+log(w2);
 shift=c1*k+c2*j;M=n1*meanAt(c1)+n2*meanAt(c2);V=n1*varianceAt(c1)+n2*varianceAt(c2);
 third=n1*thirdAt(c1)+n2*thirdAt(c2);
 for z=[0 sqrt(ell)]
  q=M-sqrt(V)*z;xx=max(q-shift,0);F=sum(exp(logW).*gammainc(xx,n,'lower'),'all');
  approximation=normcdf(-z)+third/(6*V^(3/2))*(1-z^2)*normpdf(z);
  pos=xx>0;density=zeros(size(xx));weighted=zeros(size(xx));
  density(pos)=exp(logW(pos)+(n-1)*log(xx(pos))-xx(pos)-gammaln(n));
  weighted(pos)=exp(logW(pos)+shift(pos)-q+n*log(xx(pos))-gammaln(n+1));
  f=sum(density,'all');R=sum(weighted,'all');assert(R/f<=q/n+1e-12);
  rows{end+1}=struct('n',n,'z',z,'mean_deficit',n*mu-M,'variance',V, ...
   'CDF',F,'Edgeworth_CDF',approximation,'n_scaled_CDF_error',n*(F-approximation), ...
   'weighted_low_over_density',R/f,'exact_ratio_upper_bound',q/n); %#ok<AGROW>
 end
 fprintf('n=%d two-group Edgeworth n-scaled error %.8g\n',n,rows{end}.n_scaled_CDF_error);
end
q=-1;baseCounter=integral(@(g)gammainc(q-g,2,'lower').*normpdf(g,mu,sig), ...
 -Inf,q,'AbsTol',1e-15,'RelTol',1e-12);
result=struct('matlab_version',version,'status','double diagnostics, not uniform or interval certificates', ...
 'cost_c_bounds',[0 20],'cost_tol_x',1e-10,'randomness','none', ...
 'target_cost_coefficient',-mu*(mu-1)/8,'cost_rows',{costRows}, ...
 'triangular_rows',{rows},'Gamma_base_Gaussian_splice_at_minus_one',baseCounter);
root=fileparts(fileparts(mfilename('fullpath')));fid=fopen(fullfile(root,'rounds','r022','checks.json'),'w');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
function value=meanAt(c)
value=1+c./(1+exp(c));
end
function value=varianceAt(c)
p=1./(1+exp(c));value=1+c.^2.*p.*(1-p);
end
function value=thirdAt(c)
p=1./(1+exp(c));value=2+c.^3.*p.*(1-p).*(1-2*p);
end
