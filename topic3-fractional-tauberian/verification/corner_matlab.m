started=tic; rows=[];
for t=[.02,.01,.005]
 for k=[.25,.5,1,1.5]
  h=k*t; al=1-t; mass=1-4*h*h;
  rect=[-.5+h,.5-h,-.5,.5;-.5,-.5+h,-.5+h,.5-h;.5-h,.5,-.5+h,.5-h];
  area=4*integral(@(x) arrayfun(@(z) height(z,rect,al*sqrt(mass),3*t),x),0,3*t,'AbsTol',1e-11,'RelTol',1e-7)/mass;
  rows(end+1,:)=[t,k,area,(area-2*t*t)/t^3]; %#ok<SAGROW>
 end
end
out=struct('version',version,'rows',rows,'elapsed',toc(started),'scope','adaptive original overlap numerical evidence; no certified enclosure');
fid=fopen(fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','corner_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(rows);
function y=height(x,R,threshold,bound)
lo=0;hi=bound;
for iteration=1:50
 mid=(lo+hi)/2; L=unique([2*abs(x-R(:,1));2*abs(x-R(:,2));2*abs(mid-R(:,3));2*abs(mid-R(:,4))]);L=L(L>0);v=0;
 for ell=L'
  area=sum(max(0,min(R(:,2),x+ell/2)-max(R(:,1),x-ell/2)).*max(0,min(R(:,4),mid+ell/2)-max(R(:,3),mid-ell/2)));
  v=max(v,area/ell);
 end
 if v>threshold,lo=mid;else,hi=mid;end
end
y=(lo+hi)/2;
end
