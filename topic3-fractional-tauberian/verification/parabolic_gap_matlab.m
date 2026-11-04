% Floating checks of anisotropic union and diagonal cap formulas.
started=tic; rng(20261105,'twister');count=20000;minMargin=inf;
for j=1:count
 u=.1*rand;t=u*u;a=1-t;r=t/2+t^2/2+31*t^3/32;nu=t^3/1024;eta=u^3/128;
 n=2*rand(4,1)-1;s=2*rand(4,1)-1;
 centers=[r+s(1)*eta,r-s(1)*eta;-r+s(2)*eta,r+s(2)*eta;-r-s(3)*eta,-r+s(3)*eta;r-s(4)*eta,-r-s(4)*eta]+nu*n.*[1,1;-1,1;-1,-1;1,-1];
 % Inclusion-exclusion avoids a geometric library; all 4 squares overlap.
 unionArea=0;
 for mask=1:15
  ids=find(bitget(mask,1:4));widths=a-(max(centers(ids,:),[],1)-min(centers(ids,:),[],1));
  unionArea=unionArea+(-1)^(numel(ids)+1)*prod(widths);
 end
 margin=unionArea-((a+2*r)^2-4*eta^2-16*a*nu);minMargin=min(minMargin,margin);
 assert(margin>-1e-13);
end
capRows=[];curvatureRows=[];
for t=[.01,.001,.0001,.000001]
 a=1-t;c=2*a*a-1;theta=pi/4+linspace(-t/4,t/4,1001);
 rho=sqrt(c)./(4*(sin(theta).*cos(theta)).^1.5);curvatureRows(end+1,:)=[t,min(rho),max(rho)]; %#ok<SAGROW>
 assert(all(rho>.5 & rho<1));
 for k=[.001,.1,1]
  d=k*t^3;L=sqrt(2*sqrt(2*c)+2*d);
  value=integral(@(z) sqrt(2)-z.^2./(sqrt(c+d*z.^2)+sqrt(c)),0,L,'AbsTol',1e-13,'RelTol',1e-12);
  capRows(end+1,:)=[t,k,value]; %#ok<SAGROW>
 end
end
root=fileparts(fileparts(mfilename('fullpath')));input=jsondecode(fileread(fullfile(root,'results','parabolic_gap_python.json')));err=0;
for j=1:size(capRows,1),err=max(err,abs(capRows(j,3)-str2double(input.cap_rows(j).cap_over_d32)));end
out=struct('version',version,'seed',20261105,'union_samples',count,'minimum_floating_union_margin',minMargin,'curvature_rows',curvatureRows,'cap_rows',capRows,'max_cap_cross_difference',err,'seconds',toc(started),'scope','floating checks; tiny union margins suffer cancellation, exact certificates use Fraction');
fid=fopen(fullfile(root,'results','parabolic_gap_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(out);
