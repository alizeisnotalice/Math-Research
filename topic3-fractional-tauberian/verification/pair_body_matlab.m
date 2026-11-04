started=tic;rows=[];intervalRows=[];maxCross=0;root=fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');input=jsondecode(fileread(fullfile(root,'pair_body_python.json')));
for alpha=[.71,.75,.9,.99,.999,.9999]
 t=1-alpha;x=t*(1+alpha);W=x/alpha;c=2*alpha^2-1;
 integralValue=alpha*W^2*integral(@(y) (1-y)./(alpha-W*y),0,1,'AbsTol',1e-15,'RelTol',1e-12);
 closed=x+c*log1p(-x/(1-x));R=.5*W^2;
 rows(end+1,:)=[alpha,integralValue,closed,abs(integralValue-closed),(integralValue-2*t*t)/t^3,R]; %#ok<SAGROW>
 if alpha<=.99
  for N=[16,64,256,1024]
   grid=linspace(0,W,N+1);values=alpha*(W-grid)./(alpha-grid);h=W/N;lo=h*(sum(values)-.5*(values(1)+values(end)));mid=((1:N)-.5)*h;hi=h*sum(alpha*(W-mid)./(alpha-mid));
   intervalRows(end+1,:)=[alpha,N,lo,hi]; %#ok<SAGROW>
   for j=1:numel(input.exact_area_enclosures)
    e=input.exact_area_enclosures(j);
    if abs(fractionValue(e.alpha)-alpha)<1e-12 && e.subintervals==N
     maxCross=max([maxCross,abs(lo-fractionValue(e.lower_rational)),abs(hi-fractionValue(e.upper_rational))]);
    end
   end
  end
 end
end
out=struct('version',version,'integral_rows',rows,'floating_enclosures',intervalRows,'max_enclosure_difference_from_rational',maxCross,'max_closed_integral_error',max(rows(:,4)),'elapsed_seconds',toc(started),'scope','independent floating full-interval integration; only Python rational bounds are certified');
fid=fopen(fullfile(root,'pair_body_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(rows);disp(maxCross);
function v=fractionValue(s)
p=split(string(s),'/');v=str2double(p(1));if numel(p)>1,v=v/str2double(p(2));end
end
