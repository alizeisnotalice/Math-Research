% Three-orbit full union with low-memory x-slab/2D difference grids.
% Parent runs serially; four locked nonzero steps, no baseline union query.
n3r11droot=fileparts(mfilename('fullpath'));n3r11dreceipt=jsondecode(fileread(fullfile(n3r11droot,'round011d_results.json')));n3r11dfam=jsondecode(fileread(fullfile(n3r11droot,'round011d_family.json')));
for n3r11dcase=1:numel(n3r11dreceipt.queries)
n3r11din=n3r11dreceipt.queries(n3r11dcase).parameters;n3r11din.keys=n3r11dfam.keys;n3r11dp=zeros(1,3);n3r11dph=zeros(3,3);n3r11dw=zeros(1,3);
for n3r11dj=1:3,n3r11dp(n3r11dj)=n3r11dparse(n3r11din.periods{n3r11dj});n3r11dw(n3r11dj)=n3r11dparse(n3r11din.weights{n3r11dj});end
for n3r11db=1:3,for n3r11dj=1:3,n3r11dph(n3r11db,n3r11dj)=n3r11dparse(n3r11din.phases{n3r11db}{n3r11dj});end,end
n3r11dkeys=n3r11din.keys;n3r11dboxes=zeros(size(n3r11dkeys,1),6);
for n3r11di=1:size(n3r11dkeys,1)
 n3r11dlo=reshape(n3r11dkeys(n3r11di,1:9),3,3)';n3r11dhi=reshape(n3r11dkeys(n3r11di,10:18),3,3)';n3r11dcnt=prod(max(0,n3r11dhi-n3r11dlo+1),2);n3r11dmass=sum(n3r11dcnt'.*n3r11dw);n3r11dr=n3r11dmass^(1/3)/2;n3r11dactive=n3r11dcnt>0;
 n3r11dlower=max((n3r11dhi(n3r11dactive,:)+n3r11dph(n3r11dactive,:)).*n3r11dp,[],1)-n3r11dr;n3r11dupper=min((n3r11dlo(n3r11dactive,:)+n3r11dph(n3r11dactive,:)).*n3r11dp,[],1)+n3r11dr;n3r11dlower=max(n3r11dlower,0);n3r11dupper=min(n3r11dupper,n3r11dp);n3r11dboxes(n3r11di,:)=[n3r11dlower(1),n3r11dupper(1),n3r11dlower(2),n3r11dupper(2),n3r11dlower(3),n3r11dupper(3)];
end
n3r11dboxes=n3r11dboxes(all(n3r11dboxes(:,2:2:6)>n3r11dboxes(:,1:2:5),2),:);n3r11dcuts=cell(1,3);n3r11dix=zeros(size(n3r11dboxes));
for n3r11dj=1:3
 n3r11dcuts{n3r11dj}=unique(n3r11dboxes(:,2*n3r11dj-1:2*n3r11dj));
 for n3r11di=1:size(n3r11dboxes,1),for n3r11dk=1:2,n3r11dix(n3r11di,2*n3r11dj-2+n3r11dk)=find(n3r11dcuts{n3r11dj}==n3r11dboxes(n3r11di,2*n3r11dj-2+n3r11dk),1);end,end
end
n3r11darea=diff(n3r11dcuts{2})*diff(n3r11dcuts{3})';n3r11dvol=0;
for n3r11dx=1:length(n3r11dcuts{1})-1
 n3r11ddelta=zeros(length(n3r11dcuts{2}),length(n3r11dcuts{3}));n3r11dselected=find(n3r11dboxes(:,1)<=n3r11dcuts{1}(n3r11dx)&n3r11dboxes(:,2)>=n3r11dcuts{1}(n3r11dx+1));
 for n3r11dq=1:numel(n3r11dselected)
  n3r11di=n3r11dselected(n3r11dq);n3r11da=n3r11dix(n3r11di,3);n3r11db=n3r11dix(n3r11di,4);n3r11dc=n3r11dix(n3r11di,5);n3r11dd=n3r11dix(n3r11di,6);n3r11ddelta(n3r11da,n3r11dc)=n3r11ddelta(n3r11da,n3r11dc)+1;n3r11ddelta(n3r11db,n3r11dc)=n3r11ddelta(n3r11db,n3r11dc)-1;n3r11ddelta(n3r11da,n3r11dd)=n3r11ddelta(n3r11da,n3r11dd)-1;n3r11ddelta(n3r11db,n3r11dd)=n3r11ddelta(n3r11db,n3r11dd)+1;
 end
 n3r11ddelta=cumsum(cumsum(n3r11ddelta,1),2);n3r11docc=n3r11ddelta(1:end-1,1:end-1)>0;n3r11dvol=n3r11dvol+diff(n3r11dcuts{1}(n3r11dx:n3r11dx+1))*sum(n3r11darea(n3r11docc));
end
fprintf('N3 round011d s=%s: boxes=%d totalmass=%.17g ratio=%.17g\n',n3r11dreceipt.queries(n3r11dcase).s,size(n3r11dboxes,1),sum(n3r11dw),n3r11dvol/sum(n3r11dw));
end
function n3r11dy=n3r11dparse(n3r11ds)
 n3r11dparts=split(string(n3r11ds),'/');n3r11dy=str2double(n3r11dparts(1));if numel(n3r11dparts)==2,n3r11dy=n3r11dy/str2double(n3r11dparts(2));end
end
