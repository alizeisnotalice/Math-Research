% Independent double-precision check of the explicit two-atom boundary family.
tauValues = [1/64, 1/128, 1/256];
result = struct('tau',cell(1,numel(tauValues)), 'tauR',[], 'delta',[]);
for j=1:numel(tauValues)
    tau=tauValues(j);
    f=@(u) ((1+2*exp(-u/tau)).^tau - (1-2*exp(-u/tau)).^tau ...
        + exp(-u) + (2-exp(-u/tau)).^tau)./(2+2*exp(-u/tau));
    lower=tau*log(3)+1e-14;
    [u,negValue]=fminbnd(@(v)-f(v),lower,1,optimset('TolX',1e-14));
    result(j).tau=tau; result(j).tauR=-negValue; result(j).delta=exp(-u/tau);
    fprintf('tau=%.12g tauR=%.16g delta=%.16g\n',tau,-negValue,result(j).delta);
end
fid=fopen(fullfile(fileparts(mfilename('fullpath')),'matlab_two_atom_output.json'),'w');
assert(fid>=0); fprintf(fid,'%s\n',jsonencode(result)); fclose(fid);
