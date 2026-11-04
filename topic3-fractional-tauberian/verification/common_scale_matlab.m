started=tic;rows=[];options=optimoptions('linprog','Display','none');
for a=[.01,.05,.1,.125]
 for ell=[.5,.9,.99,1,1.001,1.01,1.05,1.2,1.5,2]
  C=[a,0;-a,0;0,a;0,-a];R=[C(:,1)-ell/2,C(:,1)+ell/2,C(:,2)-ell/2,C(:,2)+ell/2];X=unique(R(:,1:2));Y=unique(R(:,3:4));capacity=zeros(15,1);
  for j=1:numel(X)-1
   for k=1:numel(Y)-1
    x=(X(j)+X(j+1))/2;y=(Y(k)+Y(k+1))/2;inside=x>R(:,1)&x<R(:,2)&y>R(:,3)&y<R(:,4);mask=sum(2.^((0:3)').*inside);
    if mask>0,capacity(mask)=capacity(mask)+(X(j+1)-X(j))*(Y(k+1)-Y(k));end
   end
  end
  A=zeros(4,15);for i=1:4,A(i,:)=bitget((1:15)',i);end
  [~,neg,flag]=linprog([zeros(15,1);-1],[ones(1,15),0;-A,ones(4,1)],[1;-.95*ell*ones(4,1)],[],[],[zeros(15,1);-inf],[capacity;inf],options);assert(flag>0);
  S=[ell^2+4*a*ell,ell^2,ell^2-4*a^2,(ell-2*a)^2];
  if S(4)>=1,M=4;elseif S(3)>=1,M=3+S(4);elseif S(2)>=1,M=2+S(4)+S(3);elseif S(1)>=1,M=1+S(4)+S(3)+S(2);else,M=4*ell^2;end
  rows(end+1,:)=[a,ell,-neg,M/4-.95*ell,abs(neg+M/4-.95*ell)]; %#ok<SAGROW>
 end
end
roots=[];
for alpha=[.91,.95,.99,.999,.9999]
 t=1-alpha;lo=0;hi=.125;
 for it=1:70
  a=(lo+hi)/2;s=sqrt(1+4*a*a);correction=4*a^4/(s*(1+2*a*a+s));
  if a-correction<t,lo=a;else,hi=a;end
 end
 a=(lo+hi)/2;s=sqrt(1+4*a*a);fourth=4*(a/t)^4/(s*(1+2*a*a+s));roots(end+1,:)=[alpha,a,fourth,(fourth-2)/t^2]; %#ok<SAGROW>
end
out=struct('version',version,'layers',rows,'roots',roots,'max_LP_error',max(rows(:,5)),'elapsed_seconds',toc(started),'scope','independent floating overlap LP and stable scalar bisection');
fid=fopen(fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','common_scale_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(out.max_LP_error);disp(roots);
