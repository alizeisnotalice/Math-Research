% Three-orbit full union with low-memory x-slab/2D difference grids.
% Main agent runs this serially; all variables have n3r8 prefix.
n3r8root=fileparts(mfilename('fullpath'));n3r8files={'round008_certificate_input.json'};
for n3r8case=1:numel(n3r8files)
n3r8in=jsondecode(fileread(fullfile(n3r8root,n3r8files{n3r8case})));n3r8p=zeros(1,3);n3r8ph=zeros(3,3);n3r8w=zeros(1,3);
for n3r8j=1:3,n3r8p(n3r8j)=n3r8parse(n3r8in.periods{n3r8j});n3r8w(n3r8j)=n3r8parse(n3r8in.weights{n3r8j});end
for n3r8b=1:3,for n3r8j=1:3,n3r8ph(n3r8b,n3r8j)=n3r8parse(n3r8in.phases{n3r8b}{n3r8j});end,end
n3r8keys=n3r8in.keys;n3r8boxes=zeros(size(n3r8keys,1),6);
for n3r8i=1:size(n3r8keys,1)
 n3r8lo=reshape(n3r8keys(n3r8i,1:9),3,3)';n3r8hi=reshape(n3r8keys(n3r8i,10:18),3,3)';n3r8cnt=prod(max(0,n3r8hi-n3r8lo+1),2);n3r8mass=sum(n3r8cnt'.*n3r8w);n3r8r=n3r8mass^(1/3)/2;n3r8active=n3r8cnt>0;
 n3r8lower=max((n3r8hi(n3r8active,:)+n3r8ph(n3r8active,:)).*n3r8p,[],1)-n3r8r;n3r8upper=min((n3r8lo(n3r8active,:)+n3r8ph(n3r8active,:)).*n3r8p,[],1)+n3r8r;n3r8lower=max(n3r8lower,0);n3r8upper=min(n3r8upper,n3r8p);n3r8boxes(n3r8i,:)=[n3r8lower(1),n3r8upper(1),n3r8lower(2),n3r8upper(2),n3r8lower(3),n3r8upper(3)];
end
n3r8boxes=n3r8boxes(all(n3r8boxes(:,2:2:6)>n3r8boxes(:,1:2:5),2),:);n3r8cuts=cell(1,3);n3r8ix=zeros(size(n3r8boxes));
for n3r8j=1:3
 n3r8cuts{n3r8j}=unique(n3r8boxes(:,2*n3r8j-1:2*n3r8j));
 for n3r8i=1:size(n3r8boxes,1),for n3r8k=1:2,n3r8ix(n3r8i,2*n3r8j-2+n3r8k)=find(n3r8cuts{n3r8j}==n3r8boxes(n3r8i,2*n3r8j-2+n3r8k),1);end,end
end
n3r8area=diff(n3r8cuts{2})*diff(n3r8cuts{3})';n3r8vol=0;
for n3r8x=1:length(n3r8cuts{1})-1
 n3r8delta=zeros(length(n3r8cuts{2}),length(n3r8cuts{3}));n3r8selected=find(n3r8boxes(:,1)<=n3r8cuts{1}(n3r8x)&n3r8boxes(:,2)>=n3r8cuts{1}(n3r8x+1));
 for n3r8q=1:numel(n3r8selected)
  n3r8i=n3r8selected(n3r8q);n3r8a=n3r8ix(n3r8i,3);n3r8b=n3r8ix(n3r8i,4);n3r8c=n3r8ix(n3r8i,5);n3r8d=n3r8ix(n3r8i,6);n3r8delta(n3r8a,n3r8c)=n3r8delta(n3r8a,n3r8c)+1;n3r8delta(n3r8b,n3r8c)=n3r8delta(n3r8b,n3r8c)-1;n3r8delta(n3r8a,n3r8d)=n3r8delta(n3r8a,n3r8d)-1;n3r8delta(n3r8b,n3r8d)=n3r8delta(n3r8b,n3r8d)+1;
 end
 n3r8delta=cumsum(cumsum(n3r8delta,1),2);n3r8occ=n3r8delta(1:end-1,1:end-1)>0;n3r8vol=n3r8vol+diff(n3r8cuts{1}(n3r8x:n3r8x+1))*sum(n3r8area(n3r8occ));
end
fprintf('N3 round008 %s: boxes=%d totalmass=%.17g ratio=%.17g\n',n3r8files{n3r8case},size(n3r8boxes,1),sum(n3r8w),n3r8vol/sum(n3r8w));
end
function n3r8y=n3r8parse(n3r8s)
 n3r8parts=split(string(n3r8s),'/');n3r8y=str2double(n3r8parts(1));if numel(n3r8parts)==2,n3r8y=n3r8y/str2double(n3r8parts(2));end
end
