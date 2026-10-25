% Three-orbit full union with low-memory x-slab/2D difference grids.
% Main agent runs this serially; all variables have n3r6 prefix.
n3r6root=fileparts(mfilename('fullpath'));n3r6files={'round006_structural_input.json','round006_certificate_input.json'};
for n3r6case=1:numel(n3r6files)
n3r6in=jsondecode(fileread(fullfile(n3r6root,n3r6files{n3r6case})));n3r6p=zeros(1,3);n3r6ph=zeros(3,3);n3r6w=zeros(1,3);
for n3r6j=1:3,n3r6p(n3r6j)=n3r6parse(n3r6in.periods{n3r6j});n3r6w(n3r6j)=n3r6parse(n3r6in.weights{n3r6j});end
for n3r6b=1:3,for n3r6j=1:3,n3r6ph(n3r6b,n3r6j)=n3r6parse(n3r6in.phases{n3r6b}{n3r6j});end,end
n3r6keys=n3r6in.keys;n3r6boxes=zeros(size(n3r6keys,1),6);
for n3r6i=1:size(n3r6keys,1)
 n3r6lo=reshape(n3r6keys(n3r6i,1:9),3,3)';n3r6hi=reshape(n3r6keys(n3r6i,10:18),3,3)';n3r6cnt=prod(max(0,n3r6hi-n3r6lo+1),2);n3r6mass=sum(n3r6cnt'.*n3r6w);n3r6r=n3r6mass^(1/3)/2;n3r6active=n3r6cnt>0;
 n3r6lower=max((n3r6hi(n3r6active,:)+n3r6ph(n3r6active,:)).*n3r6p,[],1)-n3r6r;n3r6upper=min((n3r6lo(n3r6active,:)+n3r6ph(n3r6active,:)).*n3r6p,[],1)+n3r6r;n3r6lower=max(n3r6lower,0);n3r6upper=min(n3r6upper,n3r6p);n3r6boxes(n3r6i,:)=[n3r6lower(1),n3r6upper(1),n3r6lower(2),n3r6upper(2),n3r6lower(3),n3r6upper(3)];
end
n3r6boxes=n3r6boxes(all(n3r6boxes(:,2:2:6)>n3r6boxes(:,1:2:5),2),:);n3r6cuts=cell(1,3);n3r6ix=zeros(size(n3r6boxes));
for n3r6j=1:3
 n3r6cuts{n3r6j}=unique(n3r6boxes(:,2*n3r6j-1:2*n3r6j));
 for n3r6i=1:size(n3r6boxes,1),for n3r6k=1:2,n3r6ix(n3r6i,2*n3r6j-2+n3r6k)=find(n3r6cuts{n3r6j}==n3r6boxes(n3r6i,2*n3r6j-2+n3r6k),1);end,end
end
n3r6area=diff(n3r6cuts{2})*diff(n3r6cuts{3})';n3r6vol=0;
for n3r6x=1:length(n3r6cuts{1})-1
 n3r6delta=zeros(length(n3r6cuts{2}),length(n3r6cuts{3}));n3r6selected=find(n3r6boxes(:,1)<=n3r6cuts{1}(n3r6x)&n3r6boxes(:,2)>=n3r6cuts{1}(n3r6x+1));
 for n3r6q=1:numel(n3r6selected)
  n3r6i=n3r6selected(n3r6q);n3r6a=n3r6ix(n3r6i,3);n3r6b=n3r6ix(n3r6i,4);n3r6c=n3r6ix(n3r6i,5);n3r6d=n3r6ix(n3r6i,6);n3r6delta(n3r6a,n3r6c)=n3r6delta(n3r6a,n3r6c)+1;n3r6delta(n3r6b,n3r6c)=n3r6delta(n3r6b,n3r6c)-1;n3r6delta(n3r6a,n3r6d)=n3r6delta(n3r6a,n3r6d)-1;n3r6delta(n3r6b,n3r6d)=n3r6delta(n3r6b,n3r6d)+1;
 end
 n3r6delta=cumsum(cumsum(n3r6delta,1),2);n3r6occ=n3r6delta(1:end-1,1:end-1)>0;n3r6vol=n3r6vol+diff(n3r6cuts{1}(n3r6x:n3r6x+1))*sum(n3r6area(n3r6occ));
end
fprintf('N3 round006 %s: boxes=%d totalmass=%.17g ratio=%.17g\n',n3r6files{n3r6case},size(n3r6boxes,1),sum(n3r6w),n3r6vol/sum(n3r6w));
end
function n3r6y=n3r6parse(n3r6s)
 n3r6parts=split(string(n3r6s),'/');n3r6y=str2double(n3r6parts(1));if numel(n3r6parts)==2,n3r6y=n3r6y/str2double(n3r6parts(2));end
end
