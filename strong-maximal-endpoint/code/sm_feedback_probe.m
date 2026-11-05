function result=sm_feedback_probe(step,maxN,out)
% Newton acceleration of r011 exact cell convolution. Diagnostics only.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1, step=.0125; end
if nargin<2, maxN=2048; end
if nargin<3, out=fullfile(root,'rounds','r025','feedback_checks.json'); end
mu=fzero(@(x)x-1-exp(-x),[1 2]);sigma=sqrt(mu);
q=0:step:ceil(maxN*mu+10*sqrt(maxN)); failure=tanh(q/2);
a=-expm1(-step);b=1-a/step;decay=exp(-step); rows={};maxResidual=0;
tic;
for n=2:maxN
 forcing=b*failure(2:end)+(a-b)*failure(1:end-1);
 conv=[0 filter(1,[1 -decay],forcing)];
 lo=zeros(size(q));hi=q;c=min(mu,q);
 for it=1:12
  state=q-c; [cs,ws,slope]=evaluate(state,step,conv,failure);
  r=exp(-c);res=conv+r.*cs-(1+r).*ws;
  positive=res>0;hi(positive)=c(positive);lo(~positive)=c(~positive);
  correction=res./max((1+r).*slope,realmin);
  cand=c-correction;
  bad=~isfinite(cand)|cand<lo|cand>hi;
  cand(bad)=(lo(bad)+hi(bad))/2;
  % Exact zeros and flat tails need no change.
  cand(abs(res)<2e-16)=c(abs(res)<2e-16);
  c=cand;
 end
 state=q-c;[cs,ws,~]=evaluate(state,step,conv,failure);r=exp(-c);
 res=conv+r.*cs-(1+r).*ws;
 % Refine unresolved stationarity equations by safeguarded bisection.
 bad=abs(res)>5e-15;
 if any(bad)
  qb=q(bad);lob=lo(bad);hib=hi(bad);
  for it=1:38
   cb=(lob+hib)/2;[csb,wsb,~]=evaluate(qb-cb,step,conv,failure);
   rb=exp(-cb);rr=conv(bad)+rb.*csb-(1+rb).*wsb;
   pos=rr>0;hib(pos)=cb(pos);lob(~pos)=cb(~pos);
  end
  c(bad)=(lob+hib)/2;state=q-c;[cs,ws,~]=evaluate(state,step,conv,failure);r=exp(-c);
 end
 next=(conv+r.*cs)./(1+r);next(1)=0; next=min(max(next,0),1);
 assert(min(diff(next))>-2e-11,'Bellman monotonicity failed');
 assert(min(next)>=0&&max(next)<=1);
 maxResidual=max(maxResidual,max(abs(conv+r.*cs-(1+r).*ws).*(next<.99)));
 failure=next;
 if ismember(n,[2 4 8 16 32 64 128 256 512 1024 2048])
  z=fzero(@(zz).5*erfc(zz/sqrt(2))-sigma*zz/sqrt(n),[0 5]);
  targets=[sigma*sqrt(log(n)/n),sigma*z/sqrt(n)];
  for j=1:2
   beta=targets(j);idx=find(failure>=beta,1);
   qa=q(idx-1)+step*(beta-failure(idx-1))/(failure(idx)-failure(idx-1));
   qmu=quantile(n,mu,beta,mu);
   [copt,negq]=fminbnd(@(cc)-quantile(n,cc,beta,mu),.4,2.2,optimset('TolX',1e-8,'Display','off'));
   row=struct('n',n,'target_id',j,'beta',beta,'z_root',z,'adaptive_q',qa,'mu_q',qmu,'equal_optimized_q',-negq,'equal_optimized_c',copt,'gain_mu',qa-qmu,'gain_equal_optimized',qa+negq,'gain_opt_over_log_n',(qa+negq)/log(n));
   rows{end+1}=row; %#ok<AGROW>
   fprintf('n%d h%g target%d gainMu%.9g gainOpt%.9g normalized%.9g t%.2fs\n',n,step,j,qa-qmu,qa+negq,row.gain_opt_over_log_n,toc);
  end
 end
end
result=struct('matlab',version,'step',step,'maxN',maxN,'max_stationarity_residual',maxResidual,'status','floating state interpolation; not certified; no actual C lower','seconds',toc,'rows',{rows});
fid=fopen(out,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
function q=quantile(n,c,beta,mu)
lo=0;hi=n*mu+12*sqrt(n);
for j=1:45
 q=(lo+hi)/2;[~,~,f]=sm_tail(n,exp(-c),q);
 if f<beta,lo=q;else,hi=q;end
end
q=(lo+hi)/2;
end
function [cs,ws,slope]=evaluate(s,h,conv,w)
j=min(floor(max(s,0)/h)+1,numel(w)-1);u=max(s,0)-(j-1)*h;
slope=(w(j+1)-w(j))/h;a=-expm1(-u);
cs=exp(-u).*conv(j)+a.*w(j)+max(u-a,0).*slope;ws=w(j)+u.*slope;
end
