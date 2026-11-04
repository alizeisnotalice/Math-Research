% Independent exact-overlap implementation: operator bisection on each ray.
% Quadrature and floating arithmetic are numerical evidence, not certificates.
base=fileparts(mfilename('fullpath'));
cases=[3 .5 24;3 1 24;3 1.5 24;3 2 24;3 2.5 24; ...
    4 1 12;4 2 12;4 3 12;5 1 8;5 2.5 8;5 4 8];
ts=[.01 .001];rows=zeros(size(cases,1)*numel(ts),7);k=0;
started=tic;
for j=1:size(cases,1)
    n=cases(j,1);beta=cases(j,2);q=cases(j,3);
    [z,w]=directions(n,q);
    kappa=min([1 beta n-beta]);
    for t=ts
        lo=zeros(size(z,1),1);hi=2*t/kappa+lo;
        target=log1p(-t);
        assert(all(operatorLog(hi,z,beta)<target));
        for iteration=1:54
            mid=(lo+hi)/2;
            inside=operatorLog(mid,z,beta)>target;
            lo(inside)=mid(inside);hi(~inside)=mid(~inside);
        end
        radius=(lo+hi)/2;
        residual=max(abs(operatorLog(radius,z,beta)-target));
        value=2^n*factorial(n)/n*sum(w.*(radius/t).^n);
        assert(residual<2e-14);
        k=k+1;rows(k,:)=[n beta q t value residual numel(w)];
    end
end
result=struct('kind','independent-original-overlap-operator-bisection', ...
    'version',version,'precision','IEEE754 double','random_seed',[], ...
    'bisection_steps',54,'elapsed_seconds',toc(started), ...
    'quadrature_note','Tensor Gauss-Legendre; no certified quadrature error', ...
    'columns',{{'n','beta','Gauss_order','t','volume_over_tn','log_root_residual','directions'}}, ...
    'rows',rows);
fid=fopen(fullfile(base,'..','results','higher_dimensional_matlab.json'),'w');
assert(fid~=-1);fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('Higher-dimensional overlap checks: %d integrals, max log-root residual %.4g, %.3f sec\n', ...
    k,max(rows(:,6)),result.elapsed_seconds);

function value=operatorLog(r,z,beta)
    n=size(z,2);value=-inf(size(r));
    for sign=[-1 1]
        for i=1:n
            shift=sign*2*r.*z(:,i);
            % min{1,ell,(1+ell)/2-|x_j|}, evaluated as log1p(length-1).
            lengthShift=min(min(0,shift),shift/2-r.*z);
            candidate=(beta-n)*log1p(shift)+sum(log1p(lengthShift),2);
            value=max(value,candidate);
        end
    end
end

function [z,w]=directions(n,q)
    j=(1:q-1)';off=j./sqrt(4*j.^2-1);
    [v,d]=eig(diag(off,1)+diag(off,-1),'vector');
    [d,idx]=sort(d);nodes=(d+1)/2;weights=v(1,idx)'.^2;
    dimension=n-1;ug=cell(1,dimension);wg=cell(1,dimension);
    [ug{:}]=ndgrid(nodes);[wg{:}]=ndgrid(weights);
    count=q^dimension;z=ones(count,n);w=ones(count,1);product=ones(count,1);
    for a=1:dimension
        u=ug{a}(:);product=product.*u;z(:,a+1)=product;
        w=w.*wg{a}(:).*u.^(dimension-a);
    end
    assert(abs(sum(w)-1/factorial(n-1))<2e-14);
end
