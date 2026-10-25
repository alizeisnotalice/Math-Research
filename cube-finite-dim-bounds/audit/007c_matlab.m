% Parent-serialized floating full-union cross-check for immutable round005c.
% Optionally predefine aud7cInputPath, then eval(fileread(this_script_path)).
if ~exist('aud7cInputPath','var')
 aud7cInputPath=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lower','n3','round005c_certificate_input.json');
end
aud7cInput=jsondecode(fileread(aud7cInputPath));
aud7cP=[1.630389466062437,1.072464196067522,1.129961851558576];
aud7cPhase=[0,0,0;.613275060055078,.659523766706894,.053046037414171];
aud7cW=[1,.000731338850576];
aud7cKeys=aud7cInput.witness_keys;assert(size(aud7cKeys,1)==170);
aud7cBoxes=zeros(170,6);
for aud7cI=1:170
 aud7cLo=reshape(aud7cKeys(aud7cI,1:6),3,2)';aud7cHi=reshape(aud7cKeys(aud7cI,7:12),3,2)';
 aud7cCnt=prod(max(0,aud7cHi-aud7cLo+1),2);aud7cActive=aud7cCnt>0;
 aud7cMass=aud7cCnt'*aud7cW';aud7cRad=aud7cMass^(1/3)/2;
 aud7cLower=max(0,max((aud7cHi(aud7cActive,:)+aud7cPhase(aud7cActive,:)).*aud7cP,[],1)-aud7cRad);
 aud7cUpper=min(aud7cP,min((aud7cLo(aud7cActive,:)+aud7cPhase(aud7cActive,:)).*aud7cP,[],1)+aud7cRad);
 aud7cBoxes(aud7cI,:)=[aud7cLower(1),aud7cUpper(1),aud7cLower(2),aud7cUpper(2),aud7cLower(3),aud7cUpper(3)];
end
aud7cBoxes=aud7cBoxes(all(aud7cBoxes(:,2:2:6)>aud7cBoxes(:,1:2:5),2),:);assert(size(aud7cBoxes,1)==158);
aud7cCuts=cell(1,3);aud7cCover=cell(1,3);aud7cLen=cell(1,3);
for aud7cJ=1:3
 aud7cCuts{aud7cJ}=unique(aud7cBoxes(:,2*aud7cJ-1:2*aud7cJ));aud7cLen{aud7cJ}=diff(aud7cCuts{aud7cJ});
 aud7cCover{aud7cJ}=(aud7cCuts{aud7cJ}(1:end-1)>=aud7cBoxes(:,2*aud7cJ-1)') & (aud7cCuts{aud7cJ}(2:end)<=aud7cBoxes(:,2*aud7cJ)');
end
aud7cVol=0;
for aud7cI=1:numel(aud7cLen{1})
 for aud7cJ=1:numel(aud7cLen{2})
  aud7cActive=aud7cCover{1}(aud7cI,:)&aud7cCover{2}(aud7cJ,:);
  if any(aud7cActive)
   aud7cZ=any(aud7cCover{3}(:,aud7cActive),2);
   aud7cVol=aud7cVol+aud7cLen{1}(aud7cI)*aud7cLen{2}(aud7cJ)*sum(aud7cLen{3}(aud7cZ));
  end
 end
end
aud7cRatio=aud7cVol/sum(aud7cW);
fprintf('AUDIT007c MATLAB full union: saved_keys=%d nonempty=%d volume=%.17g mass=%.17g ratio=%.17g\n',size(aud7cKeys,1),size(aud7cBoxes,1),aud7cVol,sum(aud7cW),aud7cRatio);
