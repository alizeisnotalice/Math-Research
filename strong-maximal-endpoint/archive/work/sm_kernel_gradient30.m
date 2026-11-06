function result=sm_kernel_gradient30(out)
% Independent finite-difference checks of the exact block density derivative.
if nargin<1,out='/Users/zhengzhihao/Desktop/T/work/r030_kernel_numeric.json';end
mu=fzero(@(x)x-1-exp(-x),[1 2]);rows={};
for b=[2 4 8 16 32]
 for c=[.7 mu 1.8]
  p=1/(1+exp(c));m=1+c*p;v=1+c*c*p*(1-p);
  for z=[-.9 .13 1.1]
   t=max(.1,b*m+z*sqrt(b*v));base=den(b,c,t);exact=b*p*(base-den(b-1,c,t-c));
   hs=[1e-3 1e-4 1e-5];fd=arrayfun(@(h)(den(b,c+h,t)-den(b,c-h,t))/(2*h),hs);
   rows{end+1}=struct('b',b,'c',c,'t',t,'exact',exact,'steps',hs,'finite_differences',fd,'absolute_errors',abs(fd-exact),'relative_errors',abs(fd-exact)/max(abs(exact),1e-12)); %#ok<AGROW>
  end
 end
end
result=struct('matlab',version,'status','floating finite differences, no interval certificate','rows',{rows});
fid=fopen(out,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('cases %d max abs(h1e-4) %.3g max rel(h1e-4) %.3g\n',numel(rows),max(cellfun(@(r)r.absolute_errors(2),rows)),max(cellfun(@(r)r.relative_errors(2),rows)));
end
function f=den(b,c,t)
if t<=0,f=0;return;end
j=0:b;xx=max(t-c*j,0);w=exp(gammaln(b+1)-gammaln(j+1)-gammaln(b-j+1));
f=exp(-t-b*log1p(exp(-c))-gammaln(b))*sum(w.*xx.^(b-1).*(xx>0));
end
