function result=sm_general_halo_probe(n,b2,N,panels,mesh,out)
root='/Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal/strong-maximal-endpoint';old=jsondecode(fileread(fullfile(root,'rounds','r026','checks.json')));
base=old.rows([old.rows.n]==n & [old.rows.m]==1);c=base.c;b1=n-b2;assert(min(b1,b2)>=2);tic;
M1=setup(b1,c,N,mesh,n);M2=setup(b2,c,N,mesh,n);
lo=min(M1.lr)+min(M2.lr)-n;hi=max(M1.lr)+max(M2.lr);
assert(avg(0,0,lo)>1&&avg(0,0,hi)<1);
for it=1:42,tau=(lo+hi)/2;if avg(0,0,tau)>1,lo=tau;else,hi=tau;end;end;tau=(lo+hi)/2;
a=M2.inverse(tau-M1.lr);indicator=sum(M1.pw.*exp(tau-M1.lr-a));dmat=M1.lr+M2.lr'-tau;[hh,yy]=height(dmat);active=hh>1;lp=yy+(n-1)*log1p(yy);fac=zeros(size(dmat));fac(active)=exp(lp(active)-dmat(active)).*(-expm1(-lp(active)));scaledJ=indicator+M1.pw'*(fac*M2.pw);logJ=-tau+log(scaledJ);
q1=2*M1.mean;q2=2*M2.mean;vMin=-q1-20;vMax=q2+20;radius=min([6*sqrt(min(b1,b2)),(q1+20)/2,(q2+20)/2]);outer=panels/2;assert(mod(outer,2)==0);
if b1==b2
 vv=[linspace(0,radius,panels+1),linspace(radius,vMax,outer+1)];vv(panels+2)=[];symmetry=2;
else
 vv=[linspace(vMin,-radius,outer+1),linspace(-radius,radius,panels+1),linspace(radius,vMax,outer+1)];vv([outer+2,outer+panels+3])=[];symmetry=1;
end
curve=zeros(size(vv));margins=curve;
for j=1:numel(vv)
 v=vv(j);lo=max(v-q2,-v-q1);hi=abs(v)+n*2+10;
 assert(avg(lo+v,lo-v,tau)<=1+1e-10&&avg(hi+v,hi-v,tau)>1);
 for it=1:30,dd=(lo+hi)/2;if avg(dd+v,dd-v,tau)>1,hi=dd;else,lo=dd;end;end
 curve(j)=(lo+hi)/2;margins(j)=avg(curve(j)+v+1e-5,curve(j)-v+1e-5,tau)-1;
end
ww=exp(-2*curve);
if symmetry==2
 volume=2*(simpson(ww(1:panels+1),radius/panels)+simpson(ww(panels+1:end),(vMax-radius)/outer));
else
 volume=simpson(ww(1:outer+1),(-radius-vMin)/outer)+simpson(ww(outer+1:outer+panels+1),2*radius/panels)+simpson(ww(outer+panels+1:end),(vMax-radius)/outer);
end
result=struct('matlab',version,'status','floating quadrature/root/interpolation diagnostics; no directed continuum certificate','n',n,'blocks',[b1 b2],'N',[numel(M1.t) numel(M2.t)],'panels',panels,'mesh',mesh,'c',c,'tau',tau,'logJ',logJ,'point_log_score',-logJ,'receiver_volume',volume,'full_log_score',log(volume)-logJ,'single_block_log_score',base.log_score,'gain_over_single',log(volume)-logJ-base.log_score,'mass',[sum(M1.pw) sum(M2.pw)],'calibration_residual',avg(0,0,tau)-1,'minimum_strict_grid_margin',min(margins),'analytic_receiver_tail_bound',exp(-40),'vMin',vMin,'vMax',vMax,'radius',radius,'symmetry_factor',symmetry,'seconds',toc,'v',vv,'d',curve,'weights',ww);
fid=fopen(out,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('n%d blocks%d+%d N%d V%.9g full%.10f vsSingle%+.8g margin%.3g sec%.2f\n',n,b1,b2,N,volume,result.full_log_score,result.gain_over_single,min(margins),toc);
function v=avg(u1,u2,ta)
 e1=ellAt(M1,M1.t+u1);aa=M2.inverse(ta-e1);aa(~isfinite(e1))=Inf;prob=tail(M2,aa-u2);indicator=sum(M1.pw.*prob);
 e2=ellAt(M2,M2.t+u2);h=height(e1+e2'-ta);v=indicator+M1.pw'*(max(h-1,0)*M2.pw);
end
function [h,y]=height(dd)
 h=zeros(size(dd));h(dd>=0)=1;y=zeros(size(dd));soft=dd>log(n);ys=max(expm1((dd(soft)-log(n))/(n-1)),0);
 for kk=1:5,res=(n-2)*log1p(ys)+log(n+ys)-dd(soft);der=(n-2)./(1+ys)+1./(n+ys);ys=max(0,ys-res./der);end
 y(soft)=min(ys,log(2));h(soft)=exp(y(soft));
end
end
function M=setup(b,c,N,mesh,n)
p=1/(1+exp(c));mz=1+c*p;vz=1+c*c*p*(1-p);cut=b*mz+16*sqrt(b*vz);splits=[];if b<=4,splits=c*(1:b);end
[t,w]=quad(N,cut,splits);lr=logR(b,c,t);pw=w.*exp(lr-t);
grid=unique([logspace(-12,12,3000)';(mesh:mesh:(6*n+1000))';c*(1:b)']);lgrid=zeros(size(grid));for j=1:2000:numel(grid),ix=j:min(j+1999,numel(grid));lgrid(ix)=logR(b,c,grid(ix));end
assert(all(diff(lgrid)>0));k=0:b;lw=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1)+k*log(p)+(b-k)*log1p(-p);
M=struct('b',b,'c',c,'t',t,'pw',pw,'lr',lr,'mean',b*mz,'forward',griddedInterpolant(grid,lgrid,'linear','linear'),'inverse',griddedInterpolant(lgrid,grid,'pchip','linear'),'binomialWeights',exp(lw));
end
function ell=ellAt(M,x)
ell=-Inf(size(x));pos=x>0;ell(pos)=M.forward(x(pos));
end
function a=tail(M,q)
k=0:M.b;xx=max(q-M.c*k,0);a=gammainc(xx,M.b,'upper')*M.binomialWeights';
end
function lr=logR(b,c,t)
k=0:b;lw=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1);xx=t-c*k;ll=lw+(b-1)*log(max(xx,realmin));ll(xx<=0)=-Inf;pk=max(ll,[],2);lr=pk+log(sum(exp(ll-pk),2))-gammaln(b)-b*log1p(exp(-c));
end
function [t,w]=quad(N,cut,splits)
ii=(1:15)';bb=ii./sqrt(4*ii.^2-1);[V,D]=eig(diag(bb,1)+diag(bb,-1));[x,ix]=sort(diag(D));ww=2*V(1,ix)'.^2;panels=ceil(N/16);edges=unique([linspace(0,cut,panels+1),splits]);half=diff(edges)/2;mid=(edges(1:end-1)+edges(2:end))/2;t=x*half+mid;w=ww*half;t=t(:);w=w(:);
end
function s=simpson(y,h)
assert(mod(numel(y)-1,2)==0);s=h/3*(y(1)+y(end)+4*sum(y(2:2:end-1))+2*sum(y(3:2:end-2)));
end
