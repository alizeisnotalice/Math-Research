% Stable independent root calculation of the one-dimensional branch transition.
tValues=[1e-2 1e-3 1e-4 1e-5];
rows=zeros(numel(tValues),4);
for j=1:numel(tValues)
    t=tValues(j);
    ell=-log1p(-t);
    equation=@(beta) expm1(-ell./beta)+expm1(ell./(1-beta));
    betaStar=fzero(equation,[0.2 0.7],optimset('TolX',1e-13));
    series=0.5-t/2-t^2/4-t^3/3-3*t^4/8;
    rows(j,:)=[t betaStar series abs(betaStar-series)];
end
assert(all(rows(:,4)<1e-8));
result=struct('version',version,'kind','numerical-cross-check', ...
    'columns',{{'t','beta_star_fzero','fourth_order_series','absolute_difference'}}, ...
    'rows',rows,'scope','Candidate branch equation only; exact Tauberian formula has its stated threshold restriction.');
file=fullfile(fileparts(mfilename('fullpath')),'..','results','one_dimensional_matlab.json');
fid=fopen(file,'w');assert(fid~=-1);fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
disp(rows);
