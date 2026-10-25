% Independent floating full joint union for the rationalized N3 certificate.
n3r4root=fileparts(mfilename('fullpath'));
n3r4input=jsondecode(fileread(fullfile(n3r4root,'round004_certificate_input.json')));
n3r4p=zeros(1,3);n3r4ph=zeros(2,3);n3r4w=zeros(1,2);
for n3r4j=1:3,n3r4p(n3r4j)=n3r4parse(n3r4input.periods{n3r4j});end
for n3r4b=1:2
 n3r4w(n3r4b)=n3r4parse(n3r4input.weights{n3r4b});
 for n3r4j=1:3,n3r4ph(n3r4b,n3r4j)=n3r4parse(n3r4input.phases{n3r4b}{n3r4j});end
end
n3r4keys=n3r4input.witness_keys;n3r4boxes=zeros(size(n3r4keys,1),6);
for n3r4i=1:size(n3r4keys,1)
 n3r4lo=reshape(n3r4keys(n3r4i,1:6),3,2)';n3r4hi=reshape(n3r4keys(n3r4i,7:12),3,2)';
 n3r4counts=prod(max(0,n3r4hi-n3r4lo+1),2);n3r4mass=sum(n3r4counts'.*n3r4w);n3r4radius=n3r4mass^(1/3)/2;n3r4active=n3r4counts>0;
 n3r4lower=max((n3r4hi(n3r4active,:)+n3r4ph(n3r4active,:)).*n3r4p,[],1)-n3r4radius;
 n3r4upper=min((n3r4lo(n3r4active,:)+n3r4ph(n3r4active,:)).*n3r4p,[],1)+n3r4radius;
 n3r4lower=max(n3r4lower,0);n3r4upper=min(n3r4upper,n3r4p);
 n3r4boxes(n3r4i,:)=[n3r4lower(1),n3r4upper(1),n3r4lower(2),n3r4upper(2),n3r4lower(3),n3r4upper(3)];
end
n3r4boxes=n3r4boxes(all(n3r4boxes(:,2:2:6)>n3r4boxes(:,1:2:5),2),:);
n3r4cuts=cell(1,3);for n3r4j=1:3,n3r4cuts{n3r4j}=unique(n3r4boxes(:,2*n3r4j-1:2*n3r4j));end
n3r4delta=zeros(length(n3r4cuts{1}),length(n3r4cuts{2}),length(n3r4cuts{3}));
for n3r4i=1:size(n3r4boxes,1)
 n3r4ix=zeros(3,2);
 for n3r4j=1:3,for n3r4k=1:2,n3r4ix(n3r4j,n3r4k)=find(n3r4cuts{n3r4j}==n3r4boxes(n3r4i,2*n3r4j-2+n3r4k),1);end,end
 for n3r4e=0:7
  n3r4bits=bitget(n3r4e,1:3);n3r4idx=[n3r4ix(1,n3r4bits(1)+1),n3r4ix(2,n3r4bits(2)+1),n3r4ix(3,n3r4bits(3)+1)];
  n3r4delta(n3r4idx(1),n3r4idx(2),n3r4idx(3))=n3r4delta(n3r4idx(1),n3r4idx(2),n3r4idx(3))+(-1)^sum(n3r4bits);
 end
end
for n3r4j=1:3,n3r4delta=cumsum(n3r4delta,n3r4j);end
n3r4vol=0;
for n3r4i=1:length(n3r4cuts{1})-1
 for n3r4j=1:length(n3r4cuts{2})-1
  for n3r4k=1:length(n3r4cuts{3})-1
   if n3r4delta(n3r4i,n3r4j,n3r4k)>0,n3r4vol=n3r4vol+diff(n3r4cuts{1}(n3r4i:n3r4i+1))*diff(n3r4cuts{2}(n3r4j:n3r4j+1))*diff(n3r4cuts{3}(n3r4k:n3r4k+1));end
  end
 end
end
fprintf('N3 round004 free-satellite independent MATLAB union: boxes=%d, ratio=%.17g\n',size(n3r4boxes,1),n3r4vol/sum(n3r4w));
function n3r4y=n3r4parse(n3r4s)
 n3r4parts=split(string(n3r4s),'/');n3r4y=str2double(n3r4parts(1));if numel(n3r4parts)==2,n3r4y=n3r4y/str2double(n3r4parts(2));end
end
