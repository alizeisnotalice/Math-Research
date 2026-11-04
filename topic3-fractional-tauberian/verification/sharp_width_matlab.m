started=tic;rows=[];
for alpha=[.71,.75,.9,.99,.999]
 for epsilon=[.1,.01,.001]
  w=alpha+(1-alpha)*epsilon;
  for scale=[.1,1,10]
   W=scale/w;short=scale*w;lo=0;hi=(1/alpha-alpha)*scale;
   for it=1:60
    x=(lo+hi)/2;L=unique([abs(W-2*x),W+2*x,short]);L=L(L>0);val=0;
    for ell=L
     ix=max(0,min(W/2,x+ell/2)-max(-W/2,x-ell/2));iy=min(short,ell);val=max(val,ix*iy/ell/scale);
    end
    if val>alpha,lo=x;else,hi=x;end
   end
   actual=(lo+hi)/scale;exact=w+1/w-2*alpha;rows(end+1,:)=[alpha,epsilon,scale,actual,exact,abs(actual-exact),(1/alpha-alpha)-actual]; %#ok<SAGROW>
  end
 end
end
gridRows=[];cases=[1,.25;1,.5;1,.75;2,.3;2,1;2,1.7;3,.5;3,1.5;3,2.5;5,1;5,2.5;5,4];[u,v]=ndgrid(linspace(0,1,101));
for j=1:size(cases,1)
 n=cases(j,1);beta=cases(j,2);p=n-beta;
 for alpha=[.95,.99]
  L=alpha^(1/beta);U=alpha^(-1/p);s=L+(U-L)*u;ell=s+(U-s).*v;F=(ell+s)/2+s.^(1-n).*(1-alpha*(ell.^p+s.^p));
  endpoint=max(L^(1-n)-L,U-U^(1-n));gridRows(end+1,:)=[n,beta,alpha,max(F,[],'all'),endpoint,max(F,[],'all')-endpoint]; %#ok<SAGROW>
 end
end
out=struct('version',version,'rectangle_cases',rows,'max_original_overlap_error',max(rows(:,6)),'parameter_grid',gridRows,'max_grid_excess',max(gridRows(:,6)),'elapsed_seconds',toc(started),'scope','finite original overlap and scale-invariance checks; finite grid is not a global certificate');
fid=fopen(fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','sharp_width_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp([out.max_original_overlap_error,out.max_grid_excess]);
