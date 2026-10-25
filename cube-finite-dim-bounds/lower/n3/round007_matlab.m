% Three-orbit full union with low-memory x-slab/2D difference grids.
% Main agent runs this serially; all variables have n3r7 prefix.
n3r7root=fileparts(mfilename('fullpath'));n3r7files={'round007_certificate_input.json'};
for n3r7case=1:numel(n3r7files)
n3r7in=jsondecode(fileread(fullfile(n3r7root,n3r7files{n3r7case})));n3r7p=zeros(1,3);n3r7ph=zeros(3,3);n3r7w=zeros(1,3);
for n3r7j=1:3,n3r7p(n3r7j)=n3r7parse(n3r7in.periods{n3r7j});n3r7w(n3r7j)=n3r7parse(n3r7in.weights{n3r7j});end
for n3r7b=1:3,for n3r7j=1:3,n3r7ph(n3r7b,n3r7j)=n3r7parse(n3r7in.phases{n3r7b}{n3r7j});end,end
n3r7keys=n3r7in.keys;n3r7boxes=zeros(size(n3r7keys,1),6);
for n3r7i=1:size(n3r7keys,1)
 n3r7lo=reshape(n3r7keys(n3r7i,1:9),3,3)';n3r7hi=reshape(n3r7keys(n3r7i,10:18),3,3)';n3r7cnt=prod(max(0,n3r7hi-n3r7lo+1),2);n3r7mass=sum(n3r7cnt'.*n3r7w);n3r7r=n3r7mass^(1/3)/2;n3r7active=n3r7cnt>0;
 n3r7lower=max((n3r7hi(n3r7active,:)+n3r7ph(n3r7active,:)).*n3r7p,[],1)-n3r7r;n3r7upper=min((n3r7lo(n3r7active,:)+n3r7ph(n3r7active,:)).*n3r7p,[],1)+n3r7r;n3r7lower=max(n3r7lower,0);n3r7upper=min(n3r7upper,n3r7p);n3r7boxes(n3r7i,:)=[n3r7lower(1),n3r7upper(1),n3r7lower(2),n3r7upper(2),n3r7lower(3),n3r7upper(3)];
end
n3r7boxes=n3r7boxes(all(n3r7boxes(:,2:2:6)>n3r7boxes(:,1:2:5),2),:);n3r7cuts=cell(1,3);n3r7ix=zeros(size(n3r7boxes));
for n3r7j=1:3
 n3r7cuts{n3r7j}=unique(n3r7boxes(:,2*n3r7j-1:2*n3r7j));
 for n3r7i=1:size(n3r7boxes,1),for n3r7k=1:2,n3r7ix(n3r7i,2*n3r7j-2+n3r7k)=find(n3r7cuts{n3r7j}==n3r7boxes(n3r7i,2*n3r7j-2+n3r7k),1);end,end
end
n3r7area=diff(n3r7cuts{2})*diff(n3r7cuts{3})';n3r7vol=0;
for n3r7x=1:length(n3r7cuts{1})-1
 n3r7delta=zeros(length(n3r7cuts{2}),length(n3r7cuts{3}));n3r7selected=find(n3r7boxes(:,1)<=n3r7cuts{1}(n3r7x)&n3r7boxes(:,2)>=n3r7cuts{1}(n3r7x+1));
 for n3r7q=1:numel(n3r7selected)
  n3r7i=n3r7selected(n3r7q);n3r7a=n3r7ix(n3r7i,3);n3r7b=n3r7ix(n3r7i,4);n3r7c=n3r7ix(n3r7i,5);n3r7d=n3r7ix(n3r7i,6);n3r7delta(n3r7a,n3r7c)=n3r7delta(n3r7a,n3r7c)+1;n3r7delta(n3r7b,n3r7c)=n3r7delta(n3r7b,n3r7c)-1;n3r7delta(n3r7a,n3r7d)=n3r7delta(n3r7a,n3r7d)-1;n3r7delta(n3r7b,n3r7d)=n3r7delta(n3r7b,n3r7d)+1;
 end
 n3r7delta=cumsum(cumsum(n3r7delta,1),2);n3r7occ=n3r7delta(1:end-1,1:end-1)>0;n3r7vol=n3r7vol+diff(n3r7cuts{1}(n3r7x:n3r7x+1))*sum(n3r7area(n3r7occ));
end
fprintf('N3 round007 %s: boxes=%d totalmass=%.17g ratio=%.17g\n',n3r7files{n3r7case},size(n3r7boxes,1),sum(n3r7w),n3r7vol/sum(n3r7w));
end
function n3r7y=n3r7parse(n3r7s)
 n3r7parts=split(string(n3r7s),'/');n3r7y=str2double(n3r7parts(1));if numel(n3r7parts)==2,n3r7y=n3r7y/str2double(n3r7parts(2));end
end
