started=tic;rows=[];
for alpha=[.9,.99]
 for f=[.1,.5]
  ell=alpha+f*(1/alpha-alpha);
  for rr=[.95,1.05]
   a=(1-alpha)*rr;C=[a,0;-a,0;0,a;0,-a;a/2,a/2;-a/2,a/2;a/2,-a/2;-a/2,-a/2];
   R=[C(:,1)-ell/2,C(:,1)+ell/2,C(:,2)-ell/2,C(:,2)+ell/2];
   X=unique(R(:,1:2));Y=unique(R(:,3:4));capacity=zeros(255,1);
   for j=1:numel(X)-1
    for k=1:numel(Y)-1
     x=(X(j)+X(j+1))/2;y=(Y(k)+Y(k+1))/2;
     inside=x>R(:,1)&x<R(:,2)&y>R(:,3)&y<R(:,4);mask=sum(2.^((0:7)').*inside);
     if mask>0,capacity(mask)=capacity(mask)+(X(j+1)-X(j))*(Y(k+1)-Y(k));end
    end
   end
   masks=find(capacity>1e-16);c=capacity(masks);A=zeros(8,numel(c));
   for i=1:8,A(i,:)=bitget(masks,i);end
   [pr,du]=solve(A,c,alpha*ell*ones(8,1));rows(end+1,:)=[alpha,ell,a,pr,du,abs(pr-du)]; %#ok<SAGROW>
  end
 end
end
A=[1,1,0,1,0,0;1,0,1,0,1,0;0,1,1,0,0,1];[pr,du]=solve(A,[.3;.3;.3;.1;.1;.1],.64*ones(3,1));
out=struct('version',version,'rows',rows,'example',[pr,du],'max_gap',max(rows(:,6)),'elapsed',toc(started),'scope','independent floating LP; exact certificates verified in Python');
fid=fopen(fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','witness_capacity_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(out.example);disp(out.max_gap);
function [pr,du]=solve(A,c,b)
n=numel(b);m=numel(c);opt=optimoptions('linprog','Display','none');
[~,val,flag]=linprog([zeros(m,1);-1],[ones(1,m),0;-A,ones(n,1)],[1;-b],[],[],[zeros(m,1);-inf],[c;inf],opt);assert(flag>0);pr=-val;
[~,du,flag]=linprog([-b;1;c],[A',-ones(m,1),-eye(m)],zeros(m,1),[ones(1,n),0,zeros(1,m)],1,zeros(n+1+m,1),[],opt);assert(flag>0);
end
