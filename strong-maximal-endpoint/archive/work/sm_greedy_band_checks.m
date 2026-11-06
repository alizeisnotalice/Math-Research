function result=sm_greedy_band_checks(out)
mu=fzero(@(x)x-1-exp(-x),[1 2]);p=1/(1+exp(mu));rows={};
for k=[32 64 128 256 1024 4096 16384 65536]
 L=log(k+1);sd=sqrt(k*p*(1-p));idx=max(0,floor(k*p-12*sd)):min(k,ceil(k*p+12*sd));
 lw=gammaln(k+1)-gammaln(idx+1)-gammaln(k-idx+1)+idx*log(p)+(k-idx)*log1p(-p);w=exp(lw);w=w/sum(w);
 for zz=[1 sqrt(L)]
  x=-sqrt(mu*k)*zz;s=k*mu+x;r=s+mu;
  [f,fp,~]=density(s);g0=kernel(mu,r);D=g0-cdf(k,s);a=p*D;b=p*(2*p-1)*D+p*f;assert(b>0);
  delta=-a/b;actual=g0-kernel(mu+delta,r);guaranteed=a*a/(3*b);
  cs=linspace(mu/2,3*mu/2,33)';pcs=1./(1+exp(cs));[fd,fpd,~]=density(r-cs);diff=kernel(cs,r)-cdf(k,r-cs);third=pcs.*(1-6*pcs+6*pcs.^2).*diff+pcs.*(3*pcs-2).*fd-pcs.*fpd;M=max(abs(third));
  row=struct('k',k,'global_log',L,'in_proved_horizon_domain',k>=L^4,'standardized_distance',zz,'density',f,'density_log_derivative',fp/f,'a',a,'b_over_pf',b/(p*f),'a_over_leading',a/(p*mu/2*fp),'delta_over_x_2k',delta/(x/(2*k)),'M_grid_over_density',M/f,'M_grid_times_a_over_b2',M*abs(a)/(b*b),'actual_gain',actual,'r027_gain_bound',guaranteed,'gain_over_bound',actual/guaranteed,'gain_normalized',actual/(x*x/(k*k)*f),'asymptotic_gain_normalized',p/8,'M_grid_nodes',33);
  rows{end+1}=row; %#ok<AGROW>
  fprintf('k%d z%.4f valid%d deltaRatio%.7g gainRatio%.7g norm%.7g Mcond%.7g\n',k,zz,row.in_proved_horizon_domain,row.delta_over_x_2k,row.gain_over_bound,row.gain_normalized,row.M_grid_times_a_over_b2);
 end
end
result=struct('matlab',version,'arithmetic','double exact Gamma-binomial mixtures; binomial 12sd truncation normalized; M is sampled, NOT certified supremum','rows',{rows});fid=fopen(out,'w');assert(fid>=0);cl=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
function val=cdf(shape,q)
 xx=max(q(:)-mu*idx,0);val=gammainc(xx,shape,'lower')*w';
end
function val=kernel(cc,rr)
 pc=1./(1+exp(cc(:)));val=(1-pc)*cdf(k+1,rr)+pc.*cdf(k+1,rr-cc(:));
end
function [v,vp,vpp]=density(q)
 xx=q(:)-mu*idx;pos=xx>0;pdf=zeros(size(xx));pdf(pos)=exp((k-1)*log(xx(pos))-xx(pos)-gammaln(k));fac=zeros(size(xx));fac(pos)=(k-1)./xx(pos)-1;
 sec=zeros(size(xx));sec(pos)=fac(pos).^2-(k-1)./xx(pos).^2;v=pdf*w';vp=(pdf.*fac)*w';vpp=(pdf.*sec)*w';
end
end
