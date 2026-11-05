function result=sm_block_fisher31(order,out)
% Block-parameter Fisher/KL diagnostics. Floating quadrature, not intervals.
if nargin<1,order=16;end
if nargin<2,out=fullfile(tempdir,sprintf('block_fisher31_order%d.json',order));end
mu=fzero(@(x)x-1-exp(-x),[1 2]);p=1/(1+exp(mu));rows={};tic;
for b=[2 3 4 8 16 32 64 128 256 512 1024 2048]
 [I,mass,tail]=oneFisher(b,mu,order);
 rows{end+1}=struct('b',b,'c',mu,'fisher',I,'mass',mass,'gamma_tail_bound',tail,'ratio_to_limit',I/(p*p/2),'moment_lower_bound',p*p/(2+(3-mu)/(b*mu))); %#ok<AGROW>
 fprintf('b%d I%.12g limitRatio%.8g mass%.12g\n',b,I,I/(p*p/2),mass);
end
other={};
for c=mu+[-.1 .1]
 for b=[8 32 128 512 2048]
  [I,mass,tail]=oneFisher(b,c,order);pp=1/(1+exp(c));m=1+c*pp;v=1+c*c*pp*(1-pp);mp=pp*(m-c);
  other{end+1}=struct('b',b,'c',c,'fisher',I,'mass',mass,'gamma_tail_bound',tail,'mean_information_rate',mp*mp/v,'I_over_b',I/b); %#ok<AGROW>
 end
end
klrows={};
for b=[2 4 16 64]
 for delta=[1e-2 1e-3 1e-4]
  edges=gridEdges(b,mu,delta);[t,w]=panels(order,edges);lf=logDensity(b,mu+delta,t);la=logDensity(b,mu,t);d=lf-la;phi=d+expm1(-d);small=abs(d)<1e-3;z=d(small);phi(small)=z.^2.*(1/2+z.*(-1/6+z.*(1/24+z.*(-1/120+z/720))));
  D=sum(w.*exp(lf).*phi);assert(isfinite(D) && D>=0);idx=find(cellfun(@(r)r.b==b,rows));I=rows{idx}.fisher;
  klrows{end+1}=struct('b',b,'delta',delta,'KL',D,'KL_over_half_I_delta_squared',D/(I*delta^2/2)); %#ok<AGROW>
 end
end
result=struct('matlab',version,'status','deterministic floating Gamma-binomial quadrature, not an interval certificate','order',order,'mu',mu,'predicted_stationary_limit',p*p/2,'rows',{rows},'off_stationary',{other},'KL_rows',{klrows},'seconds',toc);
fid=fopen(out,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
function [I,mass,tail]=oneFisher(b,c,order)
edges=gridEdges(b,c,0);[t,w]=panels(order,edges);lf=logDensity(b,c,t);lm=logDensity(b-1,c,t-c);p=1/(1+exp(c));sc=-b*p*expm1(lm-lf);assert(all(isfinite(sc)));I=sum(w.*exp(lf).*sc.^2);mass=sum(w.*exp(lf));tail=1e-30;assert(isfinite(I) && I>=0 && abs(mass-1)<1e-10);
end
function edges=gridEdges(b,c,delta)
p=1/(1+exp(c));m=1+c*p;v=1+c*c*p*(1-p);T=b*max(c,c+delta)+gammaincinv(1e-30,b,'upper');edges=b*m+sqrt(b*v)*(-12:.5:12);edges=edges(edges>0 & edges<T);edges=[0 edges T];
if b<=16,edges=[edges (0:b)*c (0:b)*(c+delta)];end
edges=unique(edges);edges=edges(edges>=0 & edges<=T);
end
function lf=logDensity(b,c,t)
t=t(:);lf=-Inf(size(t));good=t>0;if ~any(good),return;end
p=1/(1+exp(c));mode=floor((b+1)*p);ww=zeros(1,b+1);ww(mode+1)=1;
right=mode+1:b;ww(mode+2:end)=cumprod((b-right+1)./right*p/(1-p));left=mode-1:-1:0;ww(1:mode)=fliplr(cumprod((left+1)./(b-left)*(1-p)/p));ww=ww/sum(ww);
x=t(good)-c*(0:b);ok=x>0;zz=(max(x,realmin)-b)/b;
if b>=32,corr=1/(12*b)-1/(360*b^3)+1/(1260*b^5)-1/(1680*b^7);else,corr=gammaln(b)-(b-.5)*log(b)+b-.5*log(2*pi);end
lm=log1p(zz)-zz;small=abs(zz)<1e-3;u=zz(small);lm(small)=u.^2.*(-1/2+u.*(1/3+u.*(-1/4+u.*(1/5+u.*(-1/6+u.*(1/7+u.*(-1/8+u.*(1/9-u/10))))))));
lp=-.5*log(2*pi*b)-corr+(b-1)*lm-zz;lp(~ok)=-Inf;lp=lp+log(ww);mx=max(lp,[],2);lf(good)=mx+log(sum(exp(lp-mx),2));
end
function [t,w]=panels(n,edges)
ii=(1:n-1)';bb=ii./sqrt(4*ii.^2-1);[V,D]=eig(diag(bb,1)+diag(bb,-1));[x,ix]=sort(diag(D));ww=2*V(1,ix)'.^2;half=diff(edges)/2;mid=(edges(1:end-1)+edges(2:end))/2;t=x*half+mid;w=ww*half;t=t(:);w=w(:);
end
