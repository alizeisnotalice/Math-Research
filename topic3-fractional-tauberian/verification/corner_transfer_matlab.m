rng(20261105);started=tic;rows=[];maxOverlapError=0;
for k=[.25,.5,.625,.75,1,1.5]
 for t=[.01,.003,.001]
  h=k*t; mass=1-4*h*h;
  R=[-.5+h,.5-h,-.5,.5;-.5,-.5+h,-.5+h,.5-h;.5-h,.5,-.5+h,.5-h];
  for x=[.05,.125,.25,.4,.5]
   lo=-50;hi=50;
   for it=1:60
    q=(lo+hi)/2;r=1+t*q;v=r*[x,1-x];z=t*v;
    scales=unique([2*abs(z(1)-R(:,1));2*abs(z(1)-R(:,2));2*abs(z(2)-R(:,3));2*abs(z(2)-R(:,4))]);M=0;
    for ell=scales'
     ar=sum(max(0,min(R(:,2),z(1)+ell/2)-max(R(:,1),z(1)-ell/2)).*max(0,min(R(:,4),z(2)+ell/2)-max(R(:,3),z(2)-ell/2)));
     a=(ell-1)/t;b=max([zeros(1,2);-a*ones(1,2);v-a/2],[],1);
     U=min(k,max(0,k+v+a/2))+min(k,max(0,k-v+a/2));
     formula=prod(1-t*b)-t*t*prod(U);maxOverlapError=max(maxOverlapError,abs(ar-formula));
     M=max(M,ar/ell/sqrt(mass));
    end
    if M>1-t,lo=q;else,hi=q;end
   end
   y=1-x;m=min(x,y);aa=[-2*m,2*m,max(-2*m,min(2*m,2*(x-k))),max(-2*m,min(2*m,2*(y-k)))];
   gv=aa.^2/4+aa/2+x*y+2*k*k-(k+max(0,k-x+aa/2)).*(k+max(0,k-y+aa/2));G=max(gv);
   rows(end+1,:)=[t,k,x,(lo+hi)/2,G,abs((lo+hi)/2-G)/t]; %#ok<SAGROW>
  end
 end
end
out=struct('version',version,'rows',rows,'max_overlap_error',maxOverlapError,'max_radial_residual_over_t',max(rows(:,6)),'elapsed_seconds',toc(started),'scope','finite floating point checks; no global interval certificate');
fid=fopen(fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','corner_transfer_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp([maxOverlapError,max(rows(:,6))]);
