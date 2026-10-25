% Independent floating full joint union for the rationalized N3 certificate.
n3r5croot=fileparts(mfilename('fullpath'));
n3r5cinput=jsondecode(fileread(fullfile(n3r5croot,'round005c_certificate_input.json')));
n3r5cp=zeros(1,3);n3r5cph=zeros(2,3);n3r5cw=zeros(1,2);
for n3r5cj=1:3,n3r5cp(n3r5cj)=n3r5cparse(n3r5cinput.periods{n3r5cj});end
for n3r5cb=1:2
 n3r5cw(n3r5cb)=n3r5cparse(n3r5cinput.weights{n3r5cb});
 for n3r5cj=1:3,n3r5cph(n3r5cb,n3r5cj)=n3r5cparse(n3r5cinput.phases{n3r5cb}{n3r5cj});end
end
n3r5ckeys=n3r5cinput.witness_keys;n3r5cboxes=zeros(size(n3r5ckeys,1),6);
for n3r5ci=1:size(n3r5ckeys,1)
 n3r5clo=reshape(n3r5ckeys(n3r5ci,1:6),3,2)';n3r5chi=reshape(n3r5ckeys(n3r5ci,7:12),3,2)';
 n3r5ccounts=prod(max(0,n3r5chi-n3r5clo+1),2);n3r5cmass=sum(n3r5ccounts'.*n3r5cw);n3r5cradius=n3r5cmass^(1/3)/2;n3r5cactive=n3r5ccounts>0;
 n3r5clower=max((n3r5chi(n3r5cactive,:)+n3r5cph(n3r5cactive,:)).*n3r5cp,[],1)-n3r5cradius;
 n3r5cupper=min((n3r5clo(n3r5cactive,:)+n3r5cph(n3r5cactive,:)).*n3r5cp,[],1)+n3r5cradius;
 n3r5clower=max(n3r5clower,0);n3r5cupper=min(n3r5cupper,n3r5cp);
 n3r5cboxes(n3r5ci,:)=[n3r5clower(1),n3r5cupper(1),n3r5clower(2),n3r5cupper(2),n3r5clower(3),n3r5cupper(3)];
end
n3r5cboxes=n3r5cboxes(all(n3r5cboxes(:,2:2:6)>n3r5cboxes(:,1:2:5),2),:);
n3r5ccuts=cell(1,3);for n3r5cj=1:3,n3r5ccuts{n3r5cj}=unique(n3r5cboxes(:,2*n3r5cj-1:2*n3r5cj));end
n3r5cdelta=zeros(length(n3r5ccuts{1}),length(n3r5ccuts{2}),length(n3r5ccuts{3}));
for n3r5ci=1:size(n3r5cboxes,1)
 n3r5cix=zeros(3,2);
 for n3r5cj=1:3,for n3r5ck=1:2,n3r5cix(n3r5cj,n3r5ck)=find(n3r5ccuts{n3r5cj}==n3r5cboxes(n3r5ci,2*n3r5cj-2+n3r5ck),1);end,end
 for n3r5ce=0:7
  n3r5cbits=bitget(n3r5ce,1:3);n3r5cidx=[n3r5cix(1,n3r5cbits(1)+1),n3r5cix(2,n3r5cbits(2)+1),n3r5cix(3,n3r5cbits(3)+1)];
  n3r5cdelta(n3r5cidx(1),n3r5cidx(2),n3r5cidx(3))=n3r5cdelta(n3r5cidx(1),n3r5cidx(2),n3r5cidx(3))+(-1)^sum(n3r5cbits);
 end
end
for n3r5cj=1:3,n3r5cdelta=cumsum(n3r5cdelta,n3r5cj);end
n3r5cvol=0;
for n3r5ci=1:length(n3r5ccuts{1})-1
 for n3r5cj=1:length(n3r5ccuts{2})-1
  for n3r5ck=1:length(n3r5ccuts{3})-1
   if n3r5cdelta(n3r5ci,n3r5cj,n3r5ck)>0,n3r5cvol=n3r5cvol+diff(n3r5ccuts{1}(n3r5ci:n3r5ci+1))*diff(n3r5ccuts{2}(n3r5cj:n3r5cj+1))*diff(n3r5ccuts{3}(n3r5ck:n3r5ck+1));end
  end
 end
end
fprintf('N3 round005c free-satellite independent MATLAB union: boxes=%d, ratio=%.17g\n',size(n3r5cboxes,1),n3r5cvol/sum(n3r5cw));
function n3r5cy=n3r5cparse(n3r5cs)
 n3r5cparts=split(string(n3r5cs),'/');n3r5cy=str2double(n3r5cparts(1));if numel(n3r5cparts)==2,n3r5cy=n3r5cy/str2double(n3r5cparts(2));end
end
