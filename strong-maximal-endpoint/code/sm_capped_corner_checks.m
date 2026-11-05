function result=sm_capped_corner_checks
% Continuous rectangle-witness envelopes; rigorous lower is the Python certificate.
f=16*ones(4);f(1,1)=4;[area,parts,pieces]=witnessVolume(f);
cost=sum(f.*(1+log(f)),'all');checks=cell(0,1);
for amp=[2 4 16]
 value=witnessVolume(amp*ones(4));formula=16*(1+4*amp*log(amp));
 checks{end+1}=struct('amplitude',amp,'numerical_witness_volume',value, ...
 'constant_rectangle_exact_volume',formula,'difference',value-formula); %#ok<AGROW>
 assert(abs(value-formula)<1e-8);
end
near=cell(0,1);
for eps0=[.1 .01 .001 .0001]
 ff=(1+eps0)*ones(4);ff(1,1)=1;vv=witnessVolume(ff);cc=sum(ff.*(1+log(ff)),'all');
 base=(1/(1+eps0)+4*log1p(eps0))/(1+log1p(eps0));
 explicit=(16+4*eps0*(1-eps0^2)*(15+3*log(4/3)))/(1+15*(1+eps0)*(1+log1p(eps0)));
 near{end+1}=struct('epsilon',eps0,'envelope_ratio',vv/cc,'explicit_strip_ratio',explicit, ...
 'constant_rectangle_ratio',base,'envelope_gain_over_epsilon',(vv/cc-base)/eps0); %#ok<AGROW>
end
result=struct('matlab_version',version,'input_amplitudes',f,'witness_volume',area, ...
 'witness_parts',parts,'cost',cost,'numerical_ratio',area/cost, ...
 'constant_amplitude_16_ratio',(1/16+4*log(16))/(1+log(16)), ...
 'envelope_piece_count',numel(pieces),'constant_checks',{checks},'near_unit_rows',{near}, ...
 'status','double-precision envelope diagnostic; rational certificate is separate');
root=fileparts(fileparts(mfilename('fullpath')));fid=fopen(fullfile(root,'rounds','r023','checks.json'),'w');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('witness=%.12g ratio=%.12g constant16=%.12g\n',area,area/cost,result.constant_amplitude_16_ratio);
end
function [area,parts,pieces]=witnessVolume(F)
assert(all(F(:)>=1)&&sum(F,'all')>16);
strips=zeros(1,4);corners=zeros(1,4);pieces=struct([]);
forms={F,flipud(F),F',flipud(F')};
for k=1:4,[strips(k),pp]=strip(forms{k});for j=1:numel(pp),pp(j).orientation=k;pp(j).kind='strip';end;pieces=[pieces pp];end %#ok<AGROW>
forms={F,flipud(F),fliplr(F),rot90(F,2)};
for k=1:4,[corners(k),pp]=cornerWitness(forms{k});for j=1:numel(pp),pp(j).orientation=k;pp(j).kind='corner';end;pieces=[pieces pp];end %#ok<AGROW>
area=16+sum(strips)+sum(corners);parts=struct('interior',16,'strips',strips,'corners',corners);
end
function [area,pieces]=strip(F)
area=0;pieces=struct([]);
for j=0:3
 fs=[0 0 0 1];
 for a=0:3
  v=sum(F(a+1:4,:),1);p=[0 cumsum(v)];vc=v(j+1);
  for l=0:j
   for r=j+1:4
   fs(end+1,:)=[0 a-4+(p(r+1)-p(l+1))/(r-l) 0 1]; %#ok<AGROW>
   end
  end
  for l=0:j,fs(end+1,:)=[vc+a-4 p(j+1)-vc*j-p(l+1)-(a-4)*l 1 -l];end %#ok<AGROW>
  for r=j+1:4,fs(end+1,:)=[-vc-a+4 p(r+1)-p(j+1)+vc*j+(a-4)*r -1 r];end %#ok<AGROW>
 end
 [value,pp]=envelope(fs,j,j+1);area=area+value;pieces=[pieces pp]; %#ok<AGROW>
end
end
function [area,pieces]=cornerWitness(F)
fs=[0 0 0 1];maxExit=4;
for a=0:3
 for b=0:3
 mass=sum(F(a+1:4,b+1:4),'all');fs(end+1,:)=[b-4 mass-(b-4)*a 1 -a]; %#ok<AGROW>
 maxExit=max(maxExit,a+mass/(4-b));
 end
end
[area,pieces]=envelope(fs,4,maxExit);
end
function [area,pieces]=envelope(fs,lo,hi)
points=[lo hi];N=size(fs,1);
for i=1:N
 for j=i+1:N
 a=fs(i,:);b=fs(j,:);p=[a(1)*b(3)-b(1)*a(3),a(1)*b(4)+a(2)*b(3)-b(1)*a(4)-b(2)*a(3),a(2)*b(4)-b(2)*a(4)];
 first=find(abs(p)>1e-12,1);if isempty(first),continue;end
 rr=roots(p(first:end));rr=real(rr(abs(imag(rr))<1e-9));points=[points rr(rr>lo+1e-12 & rr<hi-1e-12)']; %#ok<AGROW>
 end
end
points=unique(sort(points));area=0;pieces=struct([]);
for j=1:numel(points)-1
 l=points(j);r=points(j+1);if r-l<1e-12,continue;end
 mid=(l+r)/2;values=(fs(:,1)*mid+fs(:,2))./(fs(:,3)*mid+fs(:,4));[v,k]=max(values);
 if v<=0,continue;end
 a=fs(k,:);cross=a(2)*a(3)-a(1)*a(4);
 if a(3)==0,value=a(1)/(2*a(4))*(r^2-l^2)+a(2)/a(4)*(r-l);
 elseif abs(cross)<1e-12,value=a(1)/a(3)*(r-l);
 else,value=a(1)/a(3)*(r-l)+cross/a(3)^2*log((a(3)*r+a(4))/(a(3)*l+a(4)));
 end
 area=area+value;piece=struct('left',l,'right',r,'coefficients',a,'integral',value,'orientation',0,'kind','');pieces=[pieces piece]; %#ok<AGROW>
end
end
