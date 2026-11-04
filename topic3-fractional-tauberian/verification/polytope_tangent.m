% Reproducible MATLAB cross-check for polytope tangent duality.
% Uses linprog for the primal/dual LPs and convhull/polyarea for 2D geometry.

scriptDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(scriptDirectory);
resultDirectory = fullfile(projectDirectory, 'results');
if ~isfolder(resultDirectory)
    mkdir(resultDirectory);
end
outputPath = fullfile(resultDirectory, 'polytope_matlab_results.txt');
fid = fopen(outputPath, 'w');
assert(fid >= 0, 'Cannot open result file.');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'MATLAB version=%s\n', version);
fprintf(fid, 'All LP checks use double precision; reported tolerances are numerical only.\n\n');

square.A = [1; 1; 1; 1];
square.nu = [1 0; 0 1; -1 0; 0 -1];
square.h = [1/2; 1/2; 1/2; 1/2];
square.V = 1;

triangle.A = [3; 3*sqrt(2); 3];
triangle.nu = [0 -1; 1/sqrt(2) 1/sqrt(2); -1 0];
triangle.h = [1; 1/sqrt(2); 1];
triangle.V = 9/2;

for beta = [0.5 1 1.5]
    emitCase(fid, 'square[-.5,.5]^2', square, beta);
end
for beta = [0.5 1 1.5]
    emitCase(fid, 'triangle[(-1,-1),(2,-1),(-1,2)]', triangle, beta);
end

diskDefect=integral(@(theta) max(-cos(theta),0),0,2*pi,'AbsTol',1e-13,'RelTol',1e-13);
assert(abs(diskDefect-2)<1e-11);
fprintf(fid,'unitDisk beta=1: directionalDefect=%.16g\n',diskDefect);
fprintf(fid,'unitDisk beta=1: areaSOverVolume=%.16g\n',pi^2/4);
clear cleanup;
fprintf('Completed six polygon LP cases and the disk integral.\n');

function emitCase(fid, name, K, beta)
    P = tangentPolygon(K, beta);
    S = polarPolygon(P);
    areaP = polygonArea(P);
    areaS = polygonArea(S);
    directions = [cos((0:63)'*2*pi/64), sin((0:63)'*2*pi/64)];
    primalDualGap = 0;
    reflectionGap = 0;
    for j = 1:size(directions, 1)
        v = directions(j, :)';
        primal = primalValue(K, beta, v);
        dual = max(P*v);
        primalDualGap = max(primalDualGap, abs(primal-dual));
        reflectionGap = max(reflectionGap, ...
            abs(primalValue(K, 2-beta, -v) - primal));
    end
    fprintf(fid, 'case=%s, beta=%.16g\n', name, beta);
    fprintf(fid, 'volume=%.16g, targetMoment=%.16g\n', K.V, (2-beta)*K.V);
    fprintf(fid, 'PVertices=%s\n', formatPoints(P));
    fprintf(fid, 'areaP=%.16g\n', areaP);
    fprintf(fid, 'SVertices=%s\n', formatPoints(S));
    fprintf(fid, 'areaS=%.16g\n', areaS);
    fprintf(fid, 'areaSOverVolume=%.16g\n', areaS/K.V);
    fprintf(fid, 'maxPrimalDualGap64Directions=%.3e\n', primalDualGap);
    fprintf(fid, 'maxReflectionGap64Directions=%.3e\n\n', reflectionGap);
end

function P = tangentPolygon(K, beta)
    m = numel(K.A);
    target = (2-beta)*K.V;
    Q = [];
    for freeIndex = 1:m
        rest = setdiff(1:m, freeIndex, 'stable');
        for mask = 0:(2^(m-1)-1)
            bits = bitget(mask, 1:m-1)';
            q = zeros(m, 1);
            q(rest) = bits;
            known = sum(K.A(rest).*K.h(rest).*bits);
            q(freeIndex) = (target-known)/(K.A(freeIndex)*K.h(freeIndex));
            if q(freeIndex) >= -1e-12 && q(freeIndex) <= 1+1e-12
                q(freeIndex) = min(1, max(0, q(freeIndex)));
                Q(end+1, :) = q'; %#ok<AGROW>
            end
        end
    end
    Q = unique(round(Q*1e12)/1e12, 'rows', 'stable');
    raw = -(K.nu'*(K.A.*Q'))'/K.V;
    raw = unique(round(raw*1e12)/1e12, 'rows', 'stable');
    index = convhull(raw(:,1), raw(:,2));
    P = raw(index(1:end-1), :);
end

function S = polarPolygon(P)
    count = size(P, 1);
    raw = zeros(count, 2);
    for i = 1:count
        j = mod(i, count) + 1;
        raw(i, :) = [P(i, :); P(j, :)]\[1; 1];
    end
    index = convhull(raw(:,1), raw(:,2));
    S = raw(index(1:end-1), :);
end

function area = polygonArea(P)
    area = polyarea(P(:,1), P(:,2));
end

function value = primalValue(K, beta, v)
    m = numel(K.A);
    f = [2-beta; K.A/K.V];
    Aineq = [-K.h, -eye(m)];
    bineq = K.nu*v;
    options = optimoptions('linprog', 'Display', 'none', 'Algorithm', 'dual-simplex');
    [~, value, exitflag] = linprog(f, Aineq, bineq, [], [], [-Inf; zeros(m,1)], [], options);
    assert(exitflag > 0, 'linprog failed.');
end

function text = formatPoints(P)
    fragments = strings(size(P,1), 1);
    for i = 1:size(P,1)
        fragments(i) = sprintf('(%.12g,%.12g)', P(i,1), P(i,2));
    end
    text = char(strjoin(fragments, '; '));
end
