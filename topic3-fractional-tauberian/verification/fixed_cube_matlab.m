% Independent double-precision contact-scale and endpoint-integral checks.
% The analytic proof, not this numerical script, establishes the theorem.
rng(20261104, 'twister');
maxExcess = 0;
sampleCount = 0;
for n = [1 2 3 5 8]
    for beta = n * [0.05 0.25 0.5 0.75 0.95]
        for trial = 1:80
            d = abs(0.98 * rand(1,n) - 0.49);
            scales = unique([1-2*d, 1+2*d]);
            bestValue = 0;
            for ell = scales
                overlap = max(0,min(min(ones(1,n),ell), (1+ell)/2-d));
                bestValue = max(bestValue, ell^(beta-n)*prod(overlap));
            end
            probeScales = exp(log(0.005) + log(1000)*rand(1,16));
            for ell = probeScales
                overlap = max(0,min(min(ones(1,n),ell), (1+ell)/2-d));
                maxExcess = max(maxExcess,ell^(beta-n)*prod(overlap)-bestValue);
                sampleCount = sampleCount+1;
            end
        end
    end
end
assert(maxExcess < 2e-12);
% Reference values are independently integrated in the active ordered chamber.
betas = [0.2 0.49 0.5 0.51 0.8 0.99 1 1.01 1.49 1.5 1.51 1.8];
rows = zeros(numel(betas)*4,5);
k = 0;
for beta = betas
    for exponent = [8 10 12 14]
        alpha = 1-2^(-exponent);
        c = 2-beta;
        if beta < 0.5
            value = expm1(log(alpha)/beta)^2;
        elseif beta < 1
            stop = -expm1(log(alpha)/beta)/2;
            value = 8*integral(@(s) 1-2*s-alpha*(1-2*s).^(1-beta),0,stop,'AbsTol',1e-14);
        elseif beta < 1.5
            stop = expm1(-log(alpha)/c)/2;
            value = 8*integral(@(s) 1-alpha*(1+2*s).^c,0,stop,'AbsTol',1e-14);
        else
            value = expm1(-log(alpha)/c)^2;
        end
        k=k+1;
        diagonal=max(-expm1(log(alpha)/beta),expm1(-log(alpha)/c))/2;
        boundaryArea=8*integral(@(s) arrayfun(@(u) boundaryWidth(u,beta,alpha),s), ...
            0,diagonal,'AbsTol',2e-14,'RelTol',2e-10);
        rows(k,:)=[beta alpha value value/(1-alpha)^2 boundaryArea];
    end
end
result=struct('kind','numerical-cross-check','version',version,'seed',20261104, ...
    'sample_count',sampleCount,'max_sampled_excess',maxExcess, ...
    'columns',{{'beta','alpha','area','area_over_t_squared','independent_boundary_area'}},'rows',rows);
outputPath=fullfile(fileparts(mfilename('fullpath')),'..','results','fixed_cube_matlab.json');
fid=fopen(outputPath,'w');
assert(fid~=-1);
fprintf(fid,'%s\n',jsonencode(result));
fclose(fid);
fprintf('MATLAB: %d scale samples, %d area cases, excess %.3g\n',sampleCount,k,maxExcess);

function width=boundaryWidth(s,beta,alpha)
    if fourValue(s,s,beta)<=alpha
        width=0;
        return;
    end
    r=fzero(@(r) fourValue(r,s,beta)-alpha,[s 0.499999999],optimset('TolX',1e-14));
    width=max(0,r-s);
end

function value=fourValue(r,s,beta)
    value=max([(1-2*r)^beta, (1-2*s)^(beta-1)*(1-s-r), ...
        (1+2*s)^(beta-2)*(1+s-r), (1+2*r)^(beta-2)]);
end
