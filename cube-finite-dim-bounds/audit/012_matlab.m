% Parent-serialized floating full-union cross-check for immutable round006_structural.
% Optionally predefine aud12InputPath, then eval(fileread(this_script_path)).
if ~exist('aud12InputPath','var')
 aud12InputPath=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lower','n3','round006_structural_input.json');
end
aud12Input=jsondecode(fileread(aud12InputPath));
aud12P=[1.630389466062437,1.072464196067522,1.129961851558576];
aud12Phase=[0,0,0;.613275060055078,.659523766706894,.053046037414171;.386724939944922,.340476233293106,.946953962585829];
aud12W=[1,.000731338850576,.000182834712644];
aud12Keys=aud12Input.keys;assert(size(aud12Keys,1)==434);
aud12Boxes=zeros(434,6);
for aud12I=1:434
 aud12Lo=reshape(aud12Keys(aud12I,1:9),3,3)';aud12Hi=reshape(aud12Keys(aud12I,10:18),3,3)';
 aud12Cnt=prod(max(0,aud12Hi-aud12Lo+1),2);aud12Active=aud12Cnt>0;
 aud12Mass=aud12Cnt'*aud12W';aud12Rad=aud12Mass^(1/3)/2;
 aud12Lower=max(0,max((aud12Hi(aud12Active,:)+aud12Phase(aud12Active,:)).*aud12P,[],1)-aud12Rad);
 aud12Upper=min(aud12P,min((aud12Lo(aud12Active,:)+aud12Phase(aud12Active,:)).*aud12P,[],1)+aud12Rad);
 aud12Boxes(aud12I,:)=[aud12Lower(1),aud12Upper(1),aud12Lower(2),aud12Upper(2),aud12Lower(3),aud12Upper(3)];
end
aud12Boxes=aud12Boxes(all(aud12Boxes(:,2:2:6)>aud12Boxes(:,1:2:5),2),:);assert(size(aud12Boxes,1)==422);
aud12Cuts=cell(1,3);aud12Cover=cell(1,3);aud12Len=cell(1,3);
for aud12J=1:3
 aud12Cuts{aud12J}=unique(aud12Boxes(:,2*aud12J-1:2*aud12J));aud12Len{aud12J}=diff(aud12Cuts{aud12J});
 aud12Cover{aud12J}=(aud12Cuts{aud12J}(1:end-1)>=aud12Boxes(:,2*aud12J-1)') & (aud12Cuts{aud12J}(2:end)<=aud12Boxes(:,2*aud12J)');
end
aud12Vol=0;
for aud12I=1:numel(aud12Len{1})
 for aud12J=1:numel(aud12Len{2})
  aud12Active=aud12Cover{1}(aud12I,:)&aud12Cover{2}(aud12J,:);
  if any(aud12Active)
   aud12Z=any(aud12Cover{3}(:,aud12Active),2);
   aud12Vol=aud12Vol+aud12Len{1}(aud12I)*aud12Len{2}(aud12J)*sum(aud12Len{3}(aud12Z));
  end
 end
end
aud12Ratio=aud12Vol/sum(aud12W);
fprintf('AUDIT012 MATLAB full union: saved_keys=%d nonempty=%d volume=%.17g mass=%.17g ratio=%.17g\n',size(aud12Keys,1),size(aud12Boxes,1),aud12Vol,sum(aud12W),aud12Ratio);
