% Independent canonical layer formula and stable one-variable rectangle optimization.
started=tic;maxLayerError=0;count=0;
for t=[.1,.05,.01]
 for lx=-1:.25:1
  for ly=-1:.25:1
   X=1+lx*t;Y=1+ly*t;boxes=[t-X/2,t+X/2,-X/2,X/2;-t-X/2,-t+X/2,-X/2,X/2;-Y/2,Y/2,t-Y/2,t+Y/2;-Y/2,Y/2,-t-Y/2,-t+Y/2];xx=unique(boxes(:,1:2));yy=unique(boxes(:,3:4));areas=zeros(1,4);
   for j=1:numel(xx)-1
    for k=1:numel(yy)-1
     xm=(xx(j)+xx(j+1))/2;ym=(yy(k)+yy(k+1))/2;n=sum(xm>boxes(:,1)&xm<boxes(:,2)&ym>boxes(:,3)&ym<boxes(:,4));
     if n>0,areas(1:n)=areas(1:n)+(xx(j+1)-xx(j))*(yy(k+1)-yy(k));end
    end
   end
   exact=[X*Y,X*X+Y*Y-X*Y-4*t*t,(X-2*t)*(Y-2*t)];maxLayerError=max(maxLayerError,max(abs(areas(2:4)-exact)));count=count+1;
  end
 end
end
root=fileparts(fileparts(mfilename('fullpath')));input=jsondecode(fileread(fullfile(root,'results','eight_witness_python.json')));rows=[];maxSigmaError=0;
for j=1:numel(input.rectangle_samples)
 e=input.rectangle_samples(j);alpha=str2double(e.alpha);c=2*alpha^2-1;d=(1+alpha)/(alpha+sqrt(c));us=alpha*(1+alpha)/(alpha+sqrt(c))^2;uhat=min(1,us);
 objective=@(u) max(u.*(1+alpha-alpha*u)./((alpha+(1-alpha)*u)*(1+alpha)),(d*d-u.^2)/(2*alpha*d));
 [u0,s0]=fminbnd(objective,0,1,optimset('TolX',1e-13,'Display','off'));candidates=[objective(0),s0,objective(1)];[sigma,ix]=min(candidates);uu=[0,u0,1];uopt=uu(ix);closed=objective(uhat);
 maxSigmaError=max(maxSigmaError,abs(sigma-str2double(e.sigma)));rows(end+1,:)=[alpha,uopt,sigma,closed]; %#ok<SAGROW>
end
Lrows=[];q=(4-sqrt(7))/6;p=2*q-q*q;K=.25+4*p;sig=7/20;
for t=[.001,.0001,.00001]
 a=1-t;z=sqrt(a*a-.75*(1-a*a));r0=(1-a*a)/(a+z);r=r0*(1-sig*t);S=sqrt(1+K*r*r);ell=S-r;
 r0MinusT=(7/4)*t*t*(1+a)/((1+z)*(a+z));eps=K*r*r/(S+1)-r0MinusT+sig*t*r0;
 Lrows(end+1,:)=[t,(ell*eps-p*r*r)/t^2,(r0MinusT-sig*t*r0)/t^2]; %#ok<SAGROW>
end
out=struct('L_corner_stable_rows',Lrows,'version',version,'exact_layer_samples',count,'max_layer_float_error',maxLayerError,'rectangle_optimization_rows',rows,'max_threshold_cross_difference',maxSigmaError,'seconds',toc(started),'scope','independent floating grid areas and bounded scalar optimization; exact mathematical statements proved separately');
fid=fopen(fullfile(root,'results','eight_witness_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(out);
