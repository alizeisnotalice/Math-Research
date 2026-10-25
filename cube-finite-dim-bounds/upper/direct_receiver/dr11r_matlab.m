% 11R: same six finite receiver configurations, no C_n claim.
% Independent vectorized breakpoint implementation; no clear/path.
upDR11Rcount=0;upDR11Rvertices=0;
for upDR11Rn=[3,4]
 for upDR11Rkind=1:3
  upDR11Rlo1=zeros(1,upDR11Rn);upDR11Rhi1=4*ones(1,upDR11Rn);
  if upDR11Rkind==1
   upDR11Rlo1(1)=-16;upDR11Rhi1(1)=0;upDR11Rrho=1/4;
  elseif upDR11Rkind==2
   upDR11Rlo1(1)=-4;upDR11Rhi1(1)=0;upDR11Rrho=1/4;
  else
   upDR11Rlo1=-4*ones(1,upDR11Rn);upDR11Rhi1=zeros(1,upDR11Rn);upDR11Rrho=4^(-upDR11Rn);
  end
  upDR11Rlo2=zeros(1,upDR11Rn);upDR11Rhi2=4*ones(1,upDR11Rn);
  upDR11Raxes=cell(1,upDR11Rn);upDR11Rmesh=cell(1,upDR11Rn);
  for upDR11Rj=1:upDR11Rn
   upDR11Raxes{upDR11Rj}=unique([upDR11Rlo1(upDR11Rj)-1,upDR11Rlo1(upDR11Rj)+1,upDR11Rhi1(upDR11Rj)-1,upDR11Rhi1(upDR11Rj)+1, ...
    upDR11Rlo2(upDR11Rj)-2,upDR11Rlo2(upDR11Rj)+2,upDR11Rhi2(upDR11Rj)-2,upDR11Rhi2(upDR11Rj)+2]);
  end
  [upDR11Rmesh{:}]=ndgrid(upDR11Raxes{:});upDR11RY=zeros(numel(upDR11Rmesh{1}),upDR11Rn);
  for upDR11Rj=1:upDR11Rn,upDR11RY(:,upDR11Rj)=upDR11Rmesh{upDR11Rj}(:);end
  upDR11Ra1=prod(max(0,min(upDR11Rhi1,upDR11RY+1)-max(upDR11Rlo1,upDR11RY-1)),2)/2^upDR11Rn;
  upDR11Ra2=prod(max(0,min(upDR11Rhi2,upDR11RY+2)-max(upDR11Rlo2,upDR11RY-2)),2)/4^upDR11Rn;
  upDR11Rh=[1-upDR11Rrho,1];upDR11Rmax=max(upDR11Rh(1)*upDR11Ra1+upDR11Ra2);
  upDR11Rrawmax=max(upDR11Ra1+upDR11Ra2);
  upDR11RV=[prod(upDR11Rhi1-upDR11Rlo1),prod(upDR11Rhi2-upDR11Rlo2)];
  upDR11Rlambda=[upDR11RV(1),upDR11RV(2)-upDR11Rrho*upDR11RV(1)];
  upDR11Rmass=upDR11RV*upDR11Rh';
  assert(upDR11Rmax==1 && all(upDR11Rlambda>=0));
  assert(upDR11Rlambda(1)==upDR11RV(1) && upDR11Rrho*upDR11Rlambda(1)+upDR11Rlambda(2)==upDR11RV(2));
  assert(sum(upDR11Rlambda)==upDR11Rmass && upDR11Rrawmax==1+upDR11Rrho);
  upDR11Rcount=upDR11Rcount+1;upDR11Rvertices=upDR11Rvertices+size(upDR11RY,1);
  fprintf('11R n%d case%d vertices%d h=(%.17g,1) maxTh %.17g rawmax %.17g mass %.17g receiver_ratio %.17g dualmass %.17g PASS\n', ...
   upDR11Rn,upDR11Rkind,size(upDR11RY,1),upDR11Rh(1),upDR11Rmax,upDR11Rrawmax,upDR11Rmass,sum(upDR11RV)/upDR11Rmass,sum(upDR11Rlambda));
 end
end
assert(upDR11Rcount==6 && upDR11Rvertices==8232);
fprintf('11R SIX_CONFIGS PASS total_vertices %d; finite receiver class only, no new C_n.\n',upDR11Rvertices);
