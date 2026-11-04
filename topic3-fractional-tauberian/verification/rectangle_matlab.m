% Independent overlap/bisection integration; finite double precision evidence.
rng(20261105); started=tic; rows=[];
for alpha=[sqrt(.5)+1e-6,.72,.8,.9,.98,.999]
 for frac=[0,.1,.5,.9,.99,.9999]
  w=1-(1-alpha)*frac; W=1/w; ss=(w+W)/2;
  A=ss-alpha; X=(1/alpha-W)/2;
  grid=unique([0,X,A]); area=0;
  for j=1:numel(grid)-1
   lo=grid(j); hi=grid(j+1);
   area=area+integral(@(x) arrayfun(@(z) height(z,w,W,alpha),x),lo,hi,'AbsTol',1e-12,'RelTol',1e-9);
  end
  area=4*area; exact=2*(ss-alpha)*(1/alpha-ss);
  rows(end+1,:)=[alpha,w,area,exact,abs(area-exact)]; %#ok<SAGROW>
 end
end
result=struct('version',version,'rows',rows,'max_absolute_error',max(rows(:,5)), 'elapsed_seconds',toc(started),'bisection_steps',60,'scope','finite floating point original overlap validation');
fid=fopen(fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','rectangle_matlab.json'),'w'); fprintf(fid,'%s',jsonencode(result)); fclose(fid); disp(result.max_absolute_error);
function y=height(x,w,W,alpha)
lo=0; hi=w/2;
for z=1:60
 mid=(lo+hi)/2;
 scales=[W-2*x,w-2*mid,W+2*x,w+2*mid]; scales=scales(scales>0);
 v=0;
 for ell=scales
  a=max(0,min(W/2,x+ell/2)-max(-W/2,x-ell/2));
  b=max(0,min(w/2,mid+ell/2)-max(-w/2,mid-ell/2));
  v=max(v,a*b/ell);
 end
 if v>alpha, lo=mid; else, hi=mid; end
end
y=(lo+hi)/2;
end
