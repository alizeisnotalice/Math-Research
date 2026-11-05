function result=sm_union_checks
% Geometry diagnostics only; exact theorem is proved in round 7.
root=fileparts(fileparts(mfilename('fullpath')));
old=rng; cleanup=onCleanup(@()rng(old)); rng(20261106,'twister');
worst=-Inf;
for j=1:10000
    count=randi([2 12]); centers=4*randn(1,count); radii=3*rand(1,count);
    b1=rand; b2=b1+2*rand;
    v1=unionLength(b1*centers,radii); v2=unionLength(b2*centers,radii);
    worst=max(worst,v1-v2);
    assert(v1<=v2+1e-12,'Uniform interval contraction failed');
end
rows=cell(0,1); q=1; shift=0.5;
for L=[4 8 12 16]
    volume=@(ell)8*exp(-ell)*(1+ell)-4*exp(-ell)*log(2/shift);
    ratio=volume(L-q)/volume(L);
    assert(ratio<=exp(q));
    rows{end+1}=struct('L',L,'q',q,'center_shift',shift, ...
        'receiver_to_input_union_ratio',ratio,'upper_bound',exp(q)); %#ok<AGROW>
    fprintf('union L=%g ratio=%.12g cap=%.12g\n',L,ratio,exp(q));
end
result=struct('matlab_version',version,'seed',20261106,'interval_cases',10000, ...
    'max_interval_contraction_excess',worst,'rows',{rows});
fid=fopen(fullfile(root,'rounds','r007','union_checks.json'),'w');
assert(fid>=0); fileCleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end

function value=unionLength(c,r)
[left,order]=sort(c-r); right=c(order)+r(order);
value=0; a=left(1); b=right(1);
for k=2:numel(c)
    if left(k)>b
        value=value+b-a; a=left(k); b=right(k);
    else
        b=max(b,right(k));
    end
end
value=value+b-a;
end
