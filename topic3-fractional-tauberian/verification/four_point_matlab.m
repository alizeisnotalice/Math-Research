% Independent floating checks; the Python rational values are the certificate.
started=tic; a=.9; r=9/160; q=.901; ep=1e-5; m=2*q-a;
g=@(z) (z-2*r).^2-2*a*z+1;
raw=[g(a),-g(q),-g(1/a),(a+2*r)^2-4*m*(m-a)-1];
robust=raw-[4*ep*(a-2*r)+4*ep^2,4*ep*(q-2*r)+4*ep^2,4*ep*(1/a-2*r)+4*ep^2,16*a*ep];
assert(all(robust>0));
root=fileparts(fileparts(mfilename('fullpath'))); data=jsondecode(fileread(fullfile(root,'results','four_point_python.json')));
maxError=0;
for j=1:numel(data.rows)
 row=data.rows(j); aa=fractionValue(row.alpha); rr=fractionValue(row.r); qq=fractionValue(row.q); ee=fractionValue(row.epsilon); mm=2*qq-aa;
 calc=[(aa-2*rr)^2-aa^2+1-aa^2-4*ee*(aa-2*rr)-4*ee^2,...
 -((qq-2*rr)^2-2*aa*qq+1)-4*ee*(qq-2*rr)-4*ee^2,...
 -((1/aa-2*rr)^2-1)-4*ee*(1/aa-2*rr)-4*ee^2,...
 (aa+2*rr)^2-4*mm*(mm-aa)-1-16*aa*ee];
 exact=[fractionValue(row.pair_margin),fractionValue(row.q_negative_margin),fractionValue(row.upper_negative_margin),fractionValue(row.mass_excess)];
 maxError=max(maxError,max(abs(calc-exact)));
end
out=struct('version',version,'raw_margins',raw,'robust_margins',robust,'checked_rows',numel(data.rows),'max_absolute_difference',maxError,'seconds',toc(started),'scope','floating cross-check; tiny margins may round to zero; exact signs come from rational Python');
fid=fopen(fullfile(root,'results','four_point_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(out);
function v=fractionValue(s)
p=split(string(s),'/');v=str2double(p(1));if numel(p)>1,v=v/str2double(p(2));end
end
