started=tic;root=fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');input=jsondecode(fileread(fullfile(root,'mesh_purification.json')));queries=input.queries_numerators_over_128/128;rows=[];maxCross=0;
for nn=1:4
 n=input.rows(nn).mesh_n;h=1/n;[ii,jj]=ndgrid(0:n-1);a=ii(:)*h;c=jj(:)*h;
 cellMass=((a+h).^2-a.^2).*((c+h).^2-c.^2)/4;right=a+cellMass/h;maxError=0;
 assert(abs(sum(cellMass)-1/4)<1e-14);
 for j=1:size(queries,1)
  x=queries(j,1);y=queries(j,2);ell=queries(j,3);
  actual=sum(max(0,min(right,x+ell)-max(a,x)).*max(0,min(c+h,y+ell)-max(c,y)));
  exact=((x+ell)^2-x^2)*((y+ell)^2-y^2)/4;err=abs(actual-exact);maxError=max(maxError,err);
  parts=split(string(input.rows(nn).errors(j)),'/');rational=str2double(parts(1));if numel(parts)>1,rational=rational/str2double(parts(2));end
  maxCross=max(maxCross,abs(err-rational));
 end
 rows(end+1,:)=[n,n*n,maxError,8*h*(1+2*h)]; %#ok<SAGROW>
end
out=struct('version',version,'rows',rows,'max_exact_python_difference',maxCross,'elapsed_seconds',toc(started),'scope','independent double geometry integrals for same 256 query parameters');
fid=fopen(fullfile(root,'mesh_purification_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(rows);disp(maxCross);
