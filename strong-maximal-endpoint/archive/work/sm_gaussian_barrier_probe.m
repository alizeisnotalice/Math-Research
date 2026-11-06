function result=sm_gaussian_barrier_probe
B=fzero(@(b)log(b)-1-1/b,[3 4]);mu=log(B);sig=sqrt(mu);C0=1e4;K=100;rows=cell(0,1);
for N=[1e4 1e6 1e8]
 L=log(N);t0=C0*L;
 for k=[1 floor(N/2) N]
  tau=k-1+t0;dh=K*L/(k+t0);b=sig*sqrt(tau);bn=sig*sqrt(tau+1);
  for z=[-6 -3 -1 0 1 3 6]*sqrt(L)
   a=mu+dh+bn*z;scale=exp(-z^2/2)/(sqrt(2*pi)*bn);
   [c,val,flag]=fminbnd(@objective,.05,5,optimset('TolX',1e-8));
   val=min([val objective(0) objective(50) objective(Inf)]);
   rows{end+1}=struct('N',N,'k',k,'z',z,'candidate_c',c,'normalized_residual',val,'exitflag',flag); %#ok<AGROW>
  end
 end
 fprintf('N=%g min normalized one-step residual %.9g\n',N,min(cellfun(@(r)r.normalized_residual,rows)));
end
result=struct('matlab_version',version,'C0',C0,'K',K,'status','sampled one-step diagnostic; not a uniform certificate','rows',{rows});
fid=fopen('/Users/zhengzhihao/Desktop/T/work/gaussian_barrier_probe.json','w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
function out=objective(c)
 if isinf(c),p=0;tail1=0;else,p=1/(1+exp(c));tail1=expcdf((a-c)/b,z>0);end
 v=(1-p)*expcdf(a/b,z>0)+p*tail1;
 if z>0,out=(normcdf(-z)-v)/scale;else,out=(v-normcdf(z))/scale;end
end
function v=expcdf(x,useUpper)
 term=.5*exp(-x.^2/2).*erfcx((b-x)/sqrt(2));
 if useUpper,v=normcdf(-x)+term;else,v=normcdf(x)-term;end
end
end
