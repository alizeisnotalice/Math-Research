% N3 round 003 independent three-orbit full joint union, floating cross-check.
n3r3root=fileparts(mfilename('fullpath'));n3r3all=jsondecode(fileread(fullfile(n3r3root,'round003_certificate_input.json')));n3r3names={'optimized','staircase','rectangular'};
for n3r3case=1:3
 n3r3in=n3r3all.(n3r3names{n3r3case});n3r3p=zeros(1,3);n3r3ph=zeros(3,3);n3r3w=zeros(1,3);
 for n3r3j=1:3,n3r3p(n3r3j)=n3r3parse(n3r3in.periods{n3r3j});n3r3w(n3r3j)=n3r3parse(n3r3in.weights{n3r3j});end
 for n3r3b=1:3,for n3r3j=1:3,n3r3ph(n3r3b,n3r3j)=n3r3parse(n3r3in.phases{n3r3b}{n3r3j});end,end
 n3r3keys=n3r3in.keys;n3r3boxes=zeros(size(n3r3keys,1),6);
 for n3r3i=1:size(n3r3keys,1)
  n3r3lo=reshape(n3r3keys(n3r3i,1:9),3,3)';n3r3hi=reshape(n3r3keys(n3r3i,10:18),3,3)';
  n3r3cnt=prod(max(0,n3r3hi-n3r3lo+1),2);n3r3mass=sum(n3r3cnt'.*n3r3w);n3r3r=n3r3mass^(1/3)/2;n3r3active=n3r3cnt>0;
  n3r3lower=max((n3r3hi(n3r3active,:)+n3r3ph(n3r3active,:)).*n3r3p,[],1)-n3r3r;n3r3upper=min((n3r3lo(n3r3active,:)+n3r3ph(n3r3active,:)).*n3r3p,[],1)+n3r3r;
  n3r3lower=max(n3r3lower,0);n3r3upper=min(n3r3upper,n3r3p);n3r3boxes(n3r3i,:)=[n3r3lower(1),n3r3upper(1),n3r3lower(2),n3r3upper(2),n3r3lower(3),n3r3upper(3)];
 end
 n3r3boxes=n3r3boxes(all(n3r3boxes(:,2:2:6)>n3r3boxes(:,1:2:5),2),:);n3r3cuts=cell(1,3);for n3r3j=1:3,n3r3cuts{n3r3j}=unique(n3r3boxes(:,2*n3r3j-1:2*n3r3j));end
 n3r3delta=zeros(length(n3r3cuts{1}),length(n3r3cuts{2}),length(n3r3cuts{3}));
 for n3r3i=1:size(n3r3boxes,1)
  n3r3ix=zeros(3,2);for n3r3j=1:3,for n3r3k=1:2,n3r3ix(n3r3j,n3r3k)=find(n3r3cuts{n3r3j}==n3r3boxes(n3r3i,2*n3r3j-2+n3r3k),1);end,end
  for n3r3e=0:7,n3r3bits=bitget(n3r3e,1:3);n3r3idx=[n3r3ix(1,n3r3bits(1)+1),n3r3ix(2,n3r3bits(2)+1),n3r3ix(3,n3r3bits(3)+1)];n3r3delta(n3r3idx(1),n3r3idx(2),n3r3idx(3))=n3r3delta(n3r3idx(1),n3r3idx(2),n3r3idx(3))+(-1)^sum(n3r3bits);end
 end
 for n3r3j=1:3,n3r3delta=cumsum(n3r3delta,n3r3j);end
 n3r3area=diff(n3r3cuts{1})*diff(n3r3cuts{2})';n3r3vol=0;for n3r3k=1:length(n3r3cuts{3})-1,n3r3occ=n3r3delta(1:end-1,1:end-1,n3r3k)>0;n3r3vol=n3r3vol+sum(n3r3area(n3r3occ))*diff(n3r3cuts{3}(n3r3k:n3r3k+1));end
 fprintf('N3 round003 %s: boxes=%d mass=%.17g ratio=%.17g\n',n3r3names{n3r3case},size(n3r3boxes,1),sum(n3r3w),n3r3vol/sum(n3r3w));
end
function n3r3y=n3r3parse(n3r3s)
 n3r3parts=split(string(n3r3s),'/');n3r3y=str2double(n3r3parts(1));if numel(n3r3parts)==2,n3r3y=n3r3y/str2double(n3r3parts(2));end
end
