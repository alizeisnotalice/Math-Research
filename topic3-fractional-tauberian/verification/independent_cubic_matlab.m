base=fullfile(fileparts(mfilename('fullpath')),'..','results');
source=jsondecode(fileread(fullfile(base,'independent_cubic_python.json')));
certMargins=zeros(1,numel(source.exact_pair_compatibility_certificates));
for j=1:numel(certMargins)
 row=source.exact_pair_compatibility_certificates(j);
 t=ratval(row.t);b=ratval(row.b);a=1-b;r=ratval(row.r);
 lo=ratval(row.low_scale);hi=ratval(row.high_scale);alpha=1-t;
 [IL,IH,bad]=pairGeometry(a,b,r,lo,hi);
 certMargins(j)=bad-(1+IL-2*alpha*lo)-(1+IH-2*alpha*hi);
 assert(certMargins(j)>0);
end
maxCornerError=0;
for j=1:numel(source.corner_regressions)
 row=source.corner_regressions(j);b=str2double(row.b);t=str2double(row.t);a=1-b;alpha=1-t;
 r=(1-alpha^2)/(alpha+sqrt(alpha^2-4*a*b*(1-alpha^2)));
 [~,~,bad]=pairGeometry(a,b,r,alpha,1/alpha);
 maxCornerError=max(maxCornerError,abs(bad/t^2-str2double(row.corner_area_over_t2)));
end
chartError=0;
for j=1:numel(source.support_chart_regressions)
 row=source.support_chart_regressions(j);t=str2double(row.t);s=str2double(row.s);a=.55;b=.45;alpha=1-t;c=2*alpha^2-1;
 r=(1-alpha^2)/(alpha+sqrt(alpha^2-4*a*b*(1-alpha^2)));p=r*[a;b];n=[alpha-2*p(2);alpha-2*p(1)];n=n/norm(n);v=[-n(2);n(1)];w=n+s*t*v;
 A=w(1)*sqrt(n(2)/n(1));B=w(2)*sqrt(n(1)/n(2));
 % Rationalized homogeneous support gap; cancellation removed before /t^2.
 f=sqrt(c)*s^2/(2*n(1)*n(2)*(sqrt(A)+sqrt(B))^2);
 chartError=max(chartError,abs(f-str2double(row.chart)));
end
assert(maxCornerError<1e-8 && chartError<1e-15);
out=struct('matlab',version,'fixed_scale_certificate_regressions',numel(certMargins), ...
 'certificate_margins',certMargins,'corner_regressions',numel(source.corner_regressions), ...
 'max_corner_coefficient_error',maxCornerError,'chart_regressions',numel(source.support_chart_regressions), ...
 'max_stable_chart_error',chartError,'scope','independent double precision geometry; universal constants supplied by analytic proof');
fid=fopen(fullfile(base,'independent_cubic_matlab.json'),'w');fprintf(fid,'%s\n',jsonencode(out));fclose(fid);disp(out);
function v=ratval(s)
 parts=split(s,'/');v=str2double(parts{1});if numel(parts)==2,v=v/str2double(parts{2});end
end
function [IL,IH,bad]=pairGeometry(a,b,r,lo,hi)
 C=r*[a,b;-a,-b;b,a;-b,-a];ell=[lo;lo;hi;hi];
 x=unique([C(:,1)-ell/2;C(:,1)+ell/2]);y=unique([C(:,2)-ell/2;C(:,2)+ell/2]);IL=0;IH=0;bad=0;
 for i=1:numel(x)-1
  for j=1:numel(y)-1
   xx=(x(i)+x(i+1))/2;yy=(y(j)+y(j+1))/2;
   inside=abs(C(:,1)-xx)<ell/2 & abs(C(:,2)-yy)<ell/2;
   area=(x(i+1)-x(i))*(y(j+1)-y(j));
   if all(inside(1:2)),IL=IL+area;end
   if all(inside(3:4)),IH=IH+area;if ~any(inside(1:2)),bad=bad+area;end;end
  end
 end
end
