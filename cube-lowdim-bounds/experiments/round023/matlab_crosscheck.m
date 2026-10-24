% Fixed-point independent MATLAB symbolic evaluation, not an interval proof.
rd = fileparts(mfilename('fullpath'));
outputPath = fullfile(rd,'matlab_result.json');
assert(~isfile(outputPath),'Existing output; do not repeat this experiment.');
src = jsondecode(fileread(fullfile(rd,'matlab_input.json')));
clockStart = tic;
records = struct([]);
costs = struct([]);
for ki = 1:numel(src.kernels)
    spec = src.kernels(ki);
    alpha = str2sym(spec.alpha); radius = str2sym(spec.R);
    step = str2sym(spec.h); singularC = str2sym(spec.c);
    expanded = zeros(0,3); % [i,j,class index], numerical indices only
    massSum = sym(0);
    for ci = 1:numel(spec.classes)
        item = spec.classes(ci); ii=item.i; jj=item.j;
        orbit = unique([ii jj;ii -jj;-ii jj;-ii -jj;jj ii;jj -ii;-jj ii;-jj -ii],'rows');
        expanded = [expanded; orbit repmat(ci,size(orbit,1),1)]; %#ok<AGROW>
        massSum = massSum + size(orbit,1)*str2sym(item.a);
        assert(toc(clockStart)<120,'Fixed calculation budget exceeded.');
    end
    cost = singularC*alpha*radius^(2-alpha)/(2-alpha)+step^2*massSum/2;
    costs(ki).name = spec.name;
    costs(ki).half_mass_80digits = char(vpa(cost,80));
    costs(ki).gap_above_3p4_80digits = char(vpa(cost-sym(17)/5,80));
    costs(ki).sum_coefficients_exact = char(massSum);
    for pi = 1:numel(src.points)
        pt = src.points(pi); uu=str2sym(pt.u); vv=str2sym(pt.v);
        assert(isAlways(uu+vv<=1));
        xx = double(uu/step)-expanded(:,1);
        yy = double(vv/step)-expanded(:,2);
        % All sample coordinates are exact dyadic rationals. This filter only
        % removes basis terms whose compact support is provably disjoint.
        active = find(abs(xx)<2 & abs(yy)<2);
        spline = sym(0);
        for ai = active'
            row=expanded(ai,:);
            bx=cubicBasis(uu/step-sym(row(1)));
            by=cubicBasis(vv/step-sym(row(2)));
            spline=spline+str2sym(spec.classes(row(3)).a)*bx*by;
        end
        value=singularC*((uu+vv)^(-alpha)-radius^(-alpha))+spline;
        rec=struct('kernel',spec.name,'point_id',pt.id,'u',pt.u,'v',pt.v,...
            'spline_exact',char(spline),'K_80digits',char(vpa(value,80)),...
            'K_minus_one_80digits',char(vpa(value-1,80)),'active_terms',numel(active));
        records(end+1)=rec; %#ok<SAGROW>
        assert(toc(clockStart)<120,'Fixed calculation budget exceeded.');
    end
end
report=struct('status','completed','K_point_evaluations',numel(records),...
 'cost_evaluations',numel(costs),'records',records,'costs',costs,...
 'elapsed_seconds',toc(clockStart),'matlab_version',version,...
 'scope','Exact rational spline plus 80-digit vpa power evaluation; numerical cross-check, not outward interval certification.');
fid=fopen(outputPath,'w'); assert(fid>=0); cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('MATLAB completed %d fixed K values; %.3f seconds.\n',numel(records),toc(clockStart));
for ki=1:numel(costs), fprintf('%s cost %s\n',costs(ki).name,costs(ki).half_mass_80digits); end
function y=cubicBasis(x)
    x=abs(x);
    if isAlways(x<=1)
        y=sym(2)/3-x^2+x^3/2;
    elseif isAlways(x<2)
        y=(2-x)^3/6;
    else
        y=sym(0);
    end
end
