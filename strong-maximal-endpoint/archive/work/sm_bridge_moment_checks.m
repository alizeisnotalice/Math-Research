function result=sm_bridge_moment_checks(order,out)
mu=fzero(@(x)x-1-exp(-x),[1 2]);sigma=sqrt(mu);rows={};cutLow=.25;cutHigh=6;C1=100;tic;
[z,qw]=zquad(order,cutLow,cutHigh);theta=cutoff(z,cutLow,cutHigh);psi=theta.*(2-theta);
for total=[1e9 1e10]
 L=log(total);c0=mu-sigma/2*sqrt(L/total);p=1/(1+exp(c0));m=1+c0*p;v=1+c0*c0*p*(1-p);
 zg=fzero(@(zz).5*erfc(zz/sqrt(2))-sqrt(v)*zz/sqrt(total),[0 sqrt(L)+1]);q=total*m-sqrt(v*total)*zg;
 eta=fzero(@(e)1/(1-e)+c0/(1+exp((1-e)*c0))-q/total,[-.01 .01]);pe=1/(1+exp((1-eta)*c0));me=1/(1-eta)+c0*pe;ve=1/(1-eta)^2+c0^2*pe*(1-pe);
 W=weights(total,p);den=density(total,c0,q,W);t0=ceil(L^4);T=floor(total/(C1*L));assert(T>t0);
 times=unique(round(exp(linspace(log(t0),log(T),3))));
 for future=times
  Wj=weights(future,p);Wp=weights(total-future,p);r=future*m+sqrt(future)*z;fj=density(future,c0,r,Wj);fp=density(total-future,c0,q-r,Wp);conditional=sqrt(future)*fj.*fp/den;
  moment=sum(qw.*z.^2.*psi.*conditional);shift=sqrt(future)*(m-me);normal=exp(-(z+shift).^2/(2*ve))/sqrt(2*pi*ve);
  row=struct('N',total,'future',future,'L',L,'c0',c0,'q',q,'eta',eta,'future_start',t0,'future_end',T,'log_window',log(T/t0),'log_window_over_L',log(T/t0)/L,'conditional_cutoff_second_moment',moment,'moment_over_mu',moment/mu,'normal_cutoff_second_moment',sum(qw.*z.^2.*psi.*normal),'mean_shift',shift,'max_relative_pdf_error_abs_z_le_3',max(abs(conditional(abs(z)<=3)./normal(abs(z)<=3)-1)),'prefix_weight_mean_error_sd',Wp.meanError,'prefix_weight_variance_ratio',Wp.varRatio,'prefix_components',numel(Wp.index),'fN_over_gaussian',den/(exp(-zg^2/2)/sqrt(2*pi*v*total)));
  rows{end+1}=row; %#ok<AGROW>
  fprintf('N%.0g j%d bridgeMoment%.10g /mu%.8g normal%.10g shift%.5g time%.1f\n',total,future,moment,moment/mu,row.normal_cutoff_second_moment,shift,toc);
 end
end
result=struct('matlab',version,'status','deterministic floating Gamma-binomial density bridge; no interval certificate','quadrature_order_per_piece',order,'z_nodes',numel(z),'cutoff',[cutLow cutHigh],'C1',C1,'Gamma_log_Stirling_remainder_bound_formula','1/(1260*shape^5)','binomial_cutoff_sd',12,'rows',{rows},'seconds',toc);
fid=fopen(out,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end
function W=weights(n,p)
sd=sqrt(n*p*(1-p));lo=max(0,floor(n*p-12*sd));hi=min(n,ceil(n*p+12*sd));mode=floor((n+1)*p);index=lo:hi;w=zeros(size(index));pivot=mode-lo+1;w(pivot)=1;
right=mode+1:hi;w(pivot+1:end)=cumprod((n-right+1)./right*p/(1-p));left=mode-1:-1:lo;wl=cumprod((left+1)./(n-left)*(1-p)/p);w(1:pivot-1)=fliplr(wl);w=w/sum(w);cent=(index-n*p)/sd;
W=struct('index',index,'w',w,'meanError',sum(w.*cent),'varRatio',sum(w.*cent.^2));
end
function values=density(shape,c,qs,W)
values=zeros(size(qs));corr=1/(12*shape)-1/(360*shape^3);base=-.5*log(2*pi*shape)-corr;
for j=1:numel(qs)
 x=qs(j)-c*W.index;good=x>0;zz=(x(good)-shape)/shape;
 lm=log1p(zz)-zz;small=abs(zz)<1e-3;u=zz(small);lm(small)=u.^2.*(-1/2+u.*(1/3+u.*(-1/4+u.*(1/5+u.*(-1/6+u.*(1/7-u/8))))));
 lp=base+(shape-1)*lm-zz;values(j)=sum(W.w(good).*exp(lp));
end
end
function th=cutoff(z,c,C)
a=abs(z);th=zeros(size(z));middle=a>=2*c&a<=C/2;th(middle)=1;up=a>c&a<2*c;down=a>C/2&a<C;
th(up)=smooth((a(up)-c)/c);th(down)=1-smooth((a(down)-C/2)/(C/2));
end
function s=smooth(u)
s=u.^4.*(35+u.*(-84+u.*(70-20*u)));
end
function [z,w]=zquad(order,c,C)
ii=(1:order-1)';bb=ii./sqrt(4*ii.^2-1);[V,D]=eig(diag(bb,1)+diag(bb,-1));[x,ix]=sort(diag(D));ww=2*V(1,ix)'.^2;edges=[-C -C/2 -2*c -c c 2*c C/2 C];z=[];w=[];
for j=[1 2 3 5 6 7]
 half=(edges(j+1)-edges(j))/2;mid=(edges(j+1)+edges(j))/2;z=[z;x*half+mid];w=[w;ww*half]; %#ok<AGROW>
end
end
