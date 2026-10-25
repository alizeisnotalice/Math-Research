% Independent floating full joint union for the rationalized N3 certificate.
n3root=fileparts(mfilename('fullpath'));
n3input=jsondecode(fileread(fullfile(n3root,'certificate_input.json')));
n3p=zeros(1,3);n3ph=zeros(2,3);n3w=zeros(1,2);
for n3j=1:3,n3p(n3j)=n3parse(n3input.periods{n3j});end
for n3b=1:2
 n3w(n3b)=n3parse(n3input.weights{n3b});
 for n3j=1:3,n3ph(n3b,n3j)=n3parse(n3input.phases{n3b}{n3j});end
end
n3keys=n3input.witness_keys;n3boxes=zeros(size(n3keys,1),6);
for n3i=1:size(n3keys,1)
 n3lo=reshape(n3keys(n3i,1:6),3,2)';n3hi=reshape(n3keys(n3i,7:12),3,2)';
 n3counts=prod(max(0,n3hi-n3lo+1),2);n3mass=sum(n3counts'.*n3w);n3radius=n3mass^(1/3)/2;n3active=n3counts>0;
 n3lower=max((n3hi(n3active,:)+n3ph(n3active,:)).*n3p,[],1)-n3radius;
 n3upper=min((n3lo(n3active,:)+n3ph(n3active,:)).*n3p,[],1)+n3radius;
 n3lower=max(n3lower,0);n3upper=min(n3upper,n3p);
 n3boxes(n3i,:)=[n3lower(1),n3upper(1),n3lower(2),n3upper(2),n3lower(3),n3upper(3)];
end
n3boxes=n3boxes(all(n3boxes(:,2:2:6)>n3boxes(:,1:2:5),2),:);
n3cuts=cell(1,3);for n3j=1:3,n3cuts{n3j}=unique(n3boxes(:,2*n3j-1:2*n3j));end
n3delta=zeros(length(n3cuts{1}),length(n3cuts{2}),length(n3cuts{3}));
for n3i=1:size(n3boxes,1)
 n3ix=zeros(3,2);
 for n3j=1:3,for n3k=1:2,n3ix(n3j,n3k)=find(n3cuts{n3j}==n3boxes(n3i,2*n3j-2+n3k),1);end,end
 for n3e=0:7
  n3bits=bitget(n3e,1:3);n3idx=[n3ix(1,n3bits(1)+1),n3ix(2,n3bits(2)+1),n3ix(3,n3bits(3)+1)];
  n3delta(n3idx(1),n3idx(2),n3idx(3))=n3delta(n3idx(1),n3idx(2),n3idx(3))+(-1)^sum(n3bits);
 end
end
for n3j=1:3,n3delta=cumsum(n3delta,n3j);end
n3vol=0;
for n3i=1:length(n3cuts{1})-1
 for n3j=1:length(n3cuts{2})-1
  for n3k=1:length(n3cuts{3})-1
   if n3delta(n3i,n3j,n3k)>0,n3vol=n3vol+diff(n3cuts{1}(n3i:n3i+1))*diff(n3cuts{2}(n3j:n3j+1))*diff(n3cuts{3}(n3k:n3k+1));end
  end
 end
end
fprintf('N3 rationalized two-orbit independent MATLAB union: boxes=%d, ratio=%.17g\n',size(n3boxes,1),n3vol/sum(n3w));
function n3y=n3parse(n3s)
 n3parts=split(string(n3s),'/');n3y=str2double(n3parts(1));if numel(n3parts)==2,n3y=n3y/str2double(n3parts(2));end
end
