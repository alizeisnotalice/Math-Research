% Independent MATLAB control for the two-atom baseline and central packets.
tauValues = [1/4, 1/2];
x = [-1, 1]; m = [4/3, 4/3];
H = subsetHalo(x, m);
pairLower = zeros(size(tauValues));
for tau = tauValues
    haloMass = weightedMass(H, tau);
    value = haloMass / sum(m .* abs(x).^(tau - 1));
    expected = 3/(4*tau)*(5/3)^tau;
    assert(abs(value-expected) < 1e-12);
    pairLower(tauValues == tau) = value;
    fprintf('tau=%g pair lower=%.16g\n', tau, value);
end
centralCosts = zeros(size(tauValues));
for tau = tauValues
    eps = 1e-6;
    centralCost = eps^(tau-1)/tau;
    centralCosts(tauValues == tau) = centralCost;
    fprintf('tau=%g central unit packet cost at eps=1e-6: %.16g\n', tau, centralCost);
end
fprintf('MATLAB %s: atomic baseline controls PASS.\n', version);
report = struct('matlab_version', version, 'tau', tauValues, ...
    'pair_lower', pairLower, 'epsilon', 1e-6, ...
    'central_unit_packet_cost', centralCosts, 'status', 'PASS');
thisFile = mfilename('fullpath');
fid = fopen(fullfile(fileparts(thisFile), 'matlab_sanity_output.json'), 'w');
assert(fid >= 0);
fprintf(fid, '%s\n', jsonencode(report));
fclose(fid);

function H = subsetHalo(x,m)
    intervals = [];
    for mask = 1:(2^numel(x)-1)
        keep = logical(bitget(mask,1:numel(x)));
        mass = sum(m(keep));
        a = max(x(keep))-mass/2; b = min(x(keep))+mass/2;
        if a < b, intervals(end+1,:)=[a,b]; end %#ok<AGROW>
    end
    intervals = sortrows(intervals,1); H=intervals(1,:);
    for j=2:size(intervals,1)
        if intervals(j,1)<=H(end,2)
            H(end,2)=max(H(end,2),intervals(j,2));
        else
            H(end+1,:)=intervals(j,:); %#ok<AGROW>
        end
    end
end

function v = weightedMass(H,tau)
    v=0;
    for j=1:size(H,1)
        a=H(j,1); b=H(j,2);
        if b<=0, v=v+((-a)^tau-(-b)^tau)/tau;
        elseif a>=0, v=v+(b^tau-a^tau)/tau;
        else, v=v+((-a)^tau+b^tau)/tau;
        end
    end
end
