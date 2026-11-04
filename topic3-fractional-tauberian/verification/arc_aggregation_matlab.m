% Independent floating arc formulas and corner-union samples.
started=tic;d=1/64;h=1/4096;k=h;rng(20261105,'twister');rows=[];randomMin=inf;
for t=[.01,.001,.0001]
 for zeta=[-.125,-.0625,-.03125,0,.03125,.0625,.125]
  a=1-t;c=2*a*a-1;theta=pi/4+zeta*t;co=cos(theta);si=sin(theta);
  rb=(a-sqrt(c*tan(theta)))/2;sb=(a-sqrt(c/tan(theta)))/2;
  r=rb-d*t^3*co;s=sb-d*t^3*si;q=a+t*t/8;M=2*q-a;mm=inf(1,4);
  for nsum=[-2,0,2]
   for ssum=[-2,0,2]
    dx=2*r+h*nsum*t^3*co-k*ssum*t^1.5*si;dy=2*s+h*nsum*t^3*si+k*ssum*t^1.5*co;
    f=@(z) (z-dx)*(z-dy)-2*a*z+1;mm(1:3)=min(mm(1:3),[f(a),-f(q),-f(1/a)]/t^3);
   end
  end
  for j=1:1000
   nn=2*rand(4,1)-1;ss=2*rand(4,1)-1;signs=[1,1;-1,1;-1,-1;1,-1];
   centers=signs.*([r,s]+h*t^3*nn.*[co,si]+k*t^1.5*ss.*[-si,co]);area=0;
   for mask=1:15
    ids=find(bitget(mask,1:4));spans=max(centers(ids,:),[],1)-min(centers(ids,:),[],1);area=area+(-1)^(numel(ids)+1)*prod(a-spans);
   end
   mm(4)=min(mm(4),(area-4*M*(M-a)-1)/t^2);
  end
  rows(end+1,:)=[t,zeta,mm]; %#ok<SAGROW>
  randomMin=min(randomMin,min(mm));
 end
end
assert(randomMin>0);root=fileparts(fileparts(mfilename('fullpath')));input=jsondecode(fileread(fullfile(root,'results','arc_aggregation_python.json')));err=0;
for j=1:size(rows,1)
 for z=1:numel(input.rows)
  rr=input.rows(z);
  if abs(str2double(rr.t)-rows(j,1))<1e-14 && abs(str2double(rr.zeta)-rows(j,2))<1e-14
   for col=1:3,err=max(err,abs(rows(j,col+2)-str2double(rr.min_normalized_margins{col})));end
  end
 end
end
out=struct('version',version,'seed',20261105,'floating_parameter_pairs',size(rows,1),'random_union_samples',21000,'rows',rows,'max_pair_cross_difference',err,'seconds',toc(started),'scope','floating diagnostics; cancellation grows at smallest t, no interval certificate');
fid=fopen(fullfile(root,'results','arc_aggregation_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(out);
