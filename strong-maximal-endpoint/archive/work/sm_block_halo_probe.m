function result=sm_block_halo_probe(n,N,vPanels,mesh,out)
% Full curved receiver diagnostics for a two-block capped likelihood profile.
root='/Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal/strong-maximal-endpoint';
old=jsondecode(fileread(fullfile(root,'rounds','r026','checks.json')));
row=old.rows([old.rows.n]==n & [old.rows.m]==2);tau0=row.tau;c=row.c;b=n/2;
p=1/(1+exp(c));mz=1+c*p;vz=1+c*c*p*(1-p);cut=b*mz+16*sqrt(b*vz);
[t,w]=quad(N,cut);lr=logR(b,c,t);pw=w.*exp(lr-t);mass=sum(pw);tic;
mainMax=6*n+1000;grid=unique([logspace(-12,12,3000)';(mesh:mesh:mainMax)';c*(1:b)']);
lgrid=zeros(size(grid));for j=1:2000:numel(grid),ii=j:min(j+1999,numel(grid));lgrid(ii)=logR(b,c,grid(ii));end
assert(all(diff(lgrid)>0));forward=griddedInterpolant(grid,lgrid,'linear','linear');inverse=griddedInterpolant(lgrid,grid,'pchip','linear');
lo=tau0-5;hi=tau0+5;assert(average(0,0,lo)>1&&average(0,0,hi)<1);
for it=1:32,tau=(lo+hi)/2;if average(0,0,tau)>1,lo=tau;else,hi=tau;end;end
tau=(lo+hi)/2;
% Cost split: exact inner reference tail plus continuous soft correction.
a=inverse(tau-lr);base=sum(pw.*exp(tau-lr-a));
d=lr+lr'-tau;[h,y]=height(d);soft=h>1;factor=zeros(size(d));lp=y+(n-1)*log1p(y);factor(soft)=exp(lp(soft)-d(soft)).*(-expm1(-lp(soft)));scaledJ=base+pw'*(factor*pw);logJ=-tau+log(scaledJ);
qUpper=2*b*mz;vMax=qUpper+20;vCentre=min(4*sqrt(n),vMax/2);outerPanels=max(8,vPanels/2);
v=[linspace(0,vCentre,vPanels+1),linspace(vCentre,vMax,outerPanels+1)];v(vPanels+2)=[];
curve=zeros(size(v));strictMargins=zeros(size(v));
for j=1:numel(v)
 vv=v(j);lo=vv-qUpper;hi=vv+n*mz+10;
 assert(average(lo+vv,lo-vv,tau)<=1+1e-10&&average(hi+vv,hi-vv,tau)>1);
 for it=1:30,dd=(lo+hi)/2;if average(dd+vv,dd-vv,tau)>1,hi=dd;else,lo=dd;end;end
 curve(j)=(lo+hi)/2;strictMargins(j)=average(curve(j)+vv+1e-5,curve(j)-vv+1e-5,tau)-1;
end
weight=exp(-2*curve);volume=2*(simpson(weight(1:vPanels+1),vCentre/vPanels)+simpson(weight(vPanels+1:end),(vMax-vCentre)/outerPanels));
result=struct('matlab',version,'status','floating quadrature/root/interpolation diagnostics; no directed average or volume certificate','n',n,'N',numel(t),'vPanels',vPanels,'mesh',mesh,'c',c,'height_cap',2,'tau_old',tau0,'tau',tau,'tau_calibration_residual',average(0,0,tau)-1,'one_block_reference_log_score',old.rows([old.rows.n]==n & [old.rows.m]==1).log_score,'point_log_score',-logJ,'receiver_volume',volume,'volume_over_sqrt_n',volume/sqrt(n),'full_log_score',log(volume)-logJ,'gain_over_single_block',log(volume)-logJ-old.rows([old.rows.n]==n & [old.rows.m]==1).log_score,'probability_mass_per_block',mass,'finite_v_max',vMax,'analytic_receiver_tail_bound',exp(-40),'minimum_grid_strict_margin',min(strictMargins),'seconds',toc,'v',v,'d',curve,'weights',weight);
fid=fopen(out,'w');assert(fid>=0);cl=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('n%d N%d V%.10g V/sqrtN%.8g point%.10f full%.10f vsSingle%+.8g tauShift%.6g margin%.3g sec%.2f\n',n,numel(t),volume,volume/sqrt(n),-logJ,result.full_log_score,result.gain_over_single_block,tau-tau0,min(strictMargins),toc);
function ell=ellAt(x)
 ell=-Inf(size(x));pos=x>0;xx=x(pos);ell(pos)=forward(xx);
end
function value=average(u1,u2,ta)
 e1=ellAt(t+u1);a2=inverse(ta-e1);a2(~isfinite(e1))=Inf;
 pp=tail(b,c,a2-u2);indicator=sum(pw.*pp);
 e2=ellAt(t+u2);dmat=e1+e2'-ta;hh=height(dmat);correction=pw'*(max(hh-1,0)*pw);value=indicator+correction;
end
function [h,y]=height(dd)
 h=zeros(size(dd));active=dd>=0;h(active)=1;y=zeros(size(dd));soft=dd>log(n);ys=max(expm1((dd(soft)-log(n))/(n-1)),0);
 for kk=1:5,rr=(n-2)*log1p(ys)+log(n+ys)-dd(soft);der=(n-2)./(1+ys)+1./(n+ys);ys=max(0,ys-rr./der);end
 y(soft)=min(ys,log(2));h(soft)=exp(y(soft));
end
end
function lr=logR(b,c,t)
k=0:b;lw=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1);xx=t-c*k;ll=lw+(b-1)*log(max(xx,realmin));ll(xx<=0)=-Inf;pk=max(ll,[],2);lr=pk+log(sum(exp(ll-pk),2))-gammaln(b)-b*log1p(exp(-c));
end
function a=tail(b,c,q)
k=0:b;p=1/(1+exp(c));lw=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1)+k*log(p)+(b-k)*log1p(-p);xx=max(q-c*k,0);a=gammainc(xx,b,'upper')*exp(lw)';
end
function [t,w]=quad(N,cut)
ii=(1:15)';bb=ii./sqrt(4*ii.^2-1);[V,D]=eig(diag(bb,1)+diag(bb,-1));[x,ix]=sort(diag(D));ww=2*V(1,ix)'.^2;panels=ceil(N/16);edges=linspace(0,cut,panels+1);half=diff(edges)/2;mid=(edges(1:end-1)+edges(2:end))/2;t=x*half+mid;w=ww*half;t=t(:);w=w(:);
end
function s=simpson(y,h)
assert(mod(numel(y)-1,2)==0);s=h/3*(y(1)+y(end)+4*sum(y(2:2:end-1))+2*sum(y(3:2:end-2)));
end
