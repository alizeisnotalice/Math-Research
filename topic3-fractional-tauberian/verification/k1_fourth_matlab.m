started=tic;rows=[];
for n=[64,128]
 beta=(1:n-1)./sqrt(4*(1:n-1).^2-1);[V,D]=eig(diag(beta,1)+diag(beta,-1));[nodes,idx]=sort(diag(D));weights=2*V(1,idx).^2;
 for t=[.02,.01,.005,.002,.001]
  R=[-.5+t,.5-t,-.5,.5;-.5,-.5+t,-.5+t,.5-t;.5-t,.5,-.5+t,.5-t];mass=1-4*t*t;intQ=0;intG=0;
  for j=1:n
   s=(nodes(j)+1)/4;lo=-40;hi=40;
   for iteration=1:60
    q=(lo+hi)/2;r=1+2*s*t+q*t*t;z=t*r*[1-s,s];
    L=unique([2*abs(z(1)-R(:,1));2*abs(z(1)-R(:,2));2*abs(z(2)-R(:,3));2*abs(z(2)-R(:,4))]);val=0;
    for ell=L'
     ar=sum(max(0,min(R(:,2),z(1)+ell/2)-max(R(:,1),z(1)-ell/2)).*max(0,min(R(:,4),z(2)+ell/2)-max(R(:,3),z(2)-ell/2)));
     val=max(val,ar/ell/sqrt(mass));
    end
    if val>1-t,lo=q;else,hi=q;end
   end
   q=(lo+hi)/2;intQ=intQ+weights(j)/4*q;intG=intG+weights(j)/4*(2*s+t*q)^2;
  end
  coef=(8*intQ+4*intG+8+8*t)/mass;rows(end+1,:)=[t,n,coef,coef+2/3]; %#ok<SAGROW>
 end
end
out=struct('version',version,'rows',rows,'elapsed_seconds',toc(started),'scope','independent original overlap double precision; unsplit Gauss not certificate');
fid=fopen(fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','k1_fourth_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(rows);
