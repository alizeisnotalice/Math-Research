% Independent Cartesian-cell layer implementation, then finite candidate oracle.
base = fullfile(fileparts(mfilename('fullpath')),'..','results');
source = jsondecode(fileread(fullfile(base,'direction_common_scale_python.json')));
maxLayerError = 0;
count = 0;
for b = (1:9)/20
    for r = [.01,.05,.2]
        for ell = [.9,1,1.1]
            got = atomLayers(ell,r,b);
            expected = closedLayers(ell,r,b);
            maxLayerError = max(maxLayerError,max(abs(got-expected)));
            count = count+1;
        end
    end
end
maxMarginError = 0;
for j = 1:numel(source.oracle_samples)
    row = source.oracle_samples(j);
    b = str2double(row.b); t = str2double(row.t);
    sigma = str2double(row.sigma); alpha = 1-t; a = 1-b;
    r0 = (1-alpha^2)/(alpha+sqrt(alpha^2-4*a*b*(1-alpha^2)));
    r = r0*(1-sigma*t);
    c = closedLayers(0,r,b);
    z = closedLayers(1,r,b)-c-1;
    candidates = [alpha,1/alpha];
    for k = 1:8
        root = (-z(k)+sqrt(z(k)^2-4*(c(k)-1)))/2;
        if root>alpha && root<1/alpha
            candidates(end+1) = root; %#ok<SAGROW>
        end
    end
    margins = zeros(size(candidates));
    for k = 1:numel(candidates)
        ell = candidates(k); v = atomLayers(ell,r,b);
        costs = (0:8)+[fliplr(cumsum(fliplr(v))),0];
        margins(k) = min(costs)-8*alpha*ell;
    end
    maxMarginError = max(maxMarginError,abs(max(margins)-str2double(row.margin)));
end
assert(maxLayerError<1e-13 && maxMarginError<1e-13);
out=struct('matlab',version,'exact_layer_float_regressions',count, ...
    'oracle_samples',numel(source.oracle_samples),'max_layer_error',maxLayerError, ...
    'max_oracle_margin_error',maxMarginError,'scope','double precision independent regressions; strict signs use exact Python intervals');
fid=fopen(fullfile(base,'direction_common_scale_matlab.json'),'w');
fprintf(fid,'%s\n',jsonencode(out));fclose(fid);disp(out);
function v=closedLayers(e,r,b)
 a=1-b;
 v=[(e+2*b*r)*(e+4*a*r-2*b*r),e^2+4*a*e*r+(-8*a*b+12*b^2)*r^2, ...
 e^2+4*b*e*r-12*b^2*r^2,e^2+4*b*e*r+(-16*a*b+4*b^2)*r^2, ...
 (e-2*a*r)*(e+2*a*r-4*b*r),(e-2*a*r)*(e+2*a*r-4*b*r),(e-2*a*r)^2,(e-2*a*r)^2];
end
function v=atomLayers(e,r,b)
 a=1-b; c=r*[a,b;-a,b;a,-b;-a,-b;b,a;-b,a;b,-a;-b,-a];
 x=unique([c(:,1)-e/2;c(:,1)+e/2]);
 y=unique([c(:,2)-e/2;c(:,2)+e/2]);v=zeros(1,8);
 for i=1:numel(x)-1
  for j=1:numel(y)-1
   xx=(x(i)+x(i+1))/2;yy=(y(j)+y(j+1))/2;
   n=sum(abs(c(:,1)-xx)<e/2 & abs(c(:,2)-yy)<e/2);
   if n>0,v(1:n)=v(1:n)+(x(i+1)-x(i))*(y(j+1)-y(j));end
  end
 end
end
