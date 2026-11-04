% Independent gamma-normalization and half-order finite-scale checks.
rng(20261104,'twister');
dimensions=[2 3 4 8 16 32 64];
logOmega=@(n) n*log(pi)/2-gammaln(1+n/2);
rows=zeros(numel(dimensions),3);
for j=1:numel(dimensions)
    n=dimensions(j);
    ellipsoid=exp(n*(logOmega(n)-logOmega(n-1)));
    cube=exp(n*log(2)-gammaln(n+1));
    assert(cube<=ellipsoid*(1+1e-12));
    rows(j,:)=[n ellipsoid cube];
end
maxGap=0;
for n=[2 3 5 8]
    for j=1:1000
        v=randn(1,n);
        candidates=[0 2*abs(v) -2*abs(v)];
        values=zeros(size(candidates));
        for k=1:numel(candidates)
            a=candidates(k);
            defect=sum(max(-v-a/2,0)+max(v-a/2,0));
            values(k)=n*a/2+defect;
        end
        maxGap=max(maxGap,abs(min(values)-sum(abs(v))));
    end
end
assert(maxGap<1e-11);
trianglePhi=@(a) a+(6*max(1-a,0)+3*max(-2-a,0))/(9/2);
assert(abs(trianglePhi(0)-4/3)<1e-14 && abs(trianglePhi(1)-1)<1e-14);
result=struct('kind','numerical-cross-check','version',version,'seed',20261104, ...
    'columns',{{'n','ellipsoid_coefficient','cube_coefficient'}},'rows',rows, ...
    'cube_directions',4000,'max_half_order_gap',maxGap, ...
    'triangle_fixed_scale',trianglePhi(0),'triangle_improved_scale',trianglePhi(1));
file=fullfile(fileparts(mfilename('fullpath')),'..','results','half_order_matlab.json');
fid=fopen(file,'w');assert(fid~=-1);fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('Half-order checks: 4000 cube directions, max gap %.3g\n',maxGap);
disp(rows);
