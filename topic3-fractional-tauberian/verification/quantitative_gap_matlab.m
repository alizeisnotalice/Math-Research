% Independent floating convex geometry and polynomial verification.
% Sampled support values are diagnostics, not certificates over all convex bodies.
started=tic; rng(20261105,'twister'); count=1500; gridCount=4096;
angle=(0:gridCount-1)'*(2*pi/gridCount); directions=[cos(angle),sin(angle)];
maxRatio=0; minSlack=inf; maxAreaIdentityError=0;
for j=1:count
    pts=randn(40,2); idx=convhull(pts); vertices=pts(idx(1:end-1),:);
    areaA=polyarea(vertices(:,1),vertices(:,2));
    diffs=reshape(permute(vertices,[1,3,2])-permute(vertices,[3,1,2]),[],2)/2;
    ii=convhull(diffs); areaS=polyarea(diffs(ii,1),diffs(ii,2)); delta=areaS-areaA;
    support=max(directions*vertices',[],2); first=2/gridCount*(support'*directions);
    odd=.5*(support-circshift(support,gridCount/2))-directions*first';
    ratio=max(abs(odd))/sqrt(delta/(2*pi)); maxRatio=max(maxRatio,ratio);minSlack=min(minSlack,1-ratio);
    % Finite Fourier energy is below the exact polygon area deficit.
    coef=fft(support)/gridCount; oddModes=3:2:255;
    fourierDelta=2*pi*sum((oddModes.^2-1).*abs(coef(oddModes+1).').^2);
    maxAreaIdentityError=max(maxAreaIdentityError,max(0,fourierDelta-delta));
end
assert(maxRatio<1+1e-5);assert(maxAreaIdentityError<1e-8);
% Six smooth perturbations verify the quadratic energy, hence sharp exponent 1/2.
smoothRows=[];
for e=[.1,.05,.01,.001,.0001,.00001]
    h=1+e*cos(3*angle);hp=-3*e*sin(3*angle);
    area=.5*(2*pi/gridCount)*sum(h.^2-hp.^2);deficit=pi-area;
    smoothRows(end+1,:)=[e,deficit,4*pi*e^2,max(abs(h-1)),deficit/e^2]; %#ok<SAGROW>
end
% Stable normalized polynomials, avoiding catastrophic cancellation of t^5 terms.
t=logspace(-12,-2,20001);
p1=(704+415*t+576*t.^2+648*t.^3)/32;
p2=(23-330*t-157*t.^2-252*t.^3-324*t.^4)/16;
p3=(24+64*t-232*t.^2+135*t.^3-66*t.^4-13*t.^5+396*t.^6-324*t.^7)/16;
p4=(304+640*t+240*t.^2+801*t.^3+1080*t.^4+1296*t.^5)/64;
mins=[min(p1),min(p2),min(p3),min(p4)];assert(all(mins>0));
momentRows=[];
for p=[64,128,256,512,1024,4096,16384]
 scaled=integral(@(u) exp((p-1)*log1p(-u/p)).*u.^6,0,p/100,'AbsTol',1e-13,'RelTol',1e-12)/768;
 fullBeta=(15/16)*prod(p./(p+(1:6)));
 momentRows(end+1,:)=[p,scaled,fullBeta]; %#ok<SAGROW>
end
root=fileparts(fileparts(mfilename('fullpath')));
out=struct('version',version,'seed',20261105,'random_convex_polygons',count,'support_grid',gridCount,'max_stability_ratio',maxRatio,'minimum_relative_slack',minSlack,'max_finite_energy_excess',maxAreaIdentityError,'smooth_perturbation_rows',smoothRows,'normalized_polynomial_samples',numel(t),'minimum_normalized_margins',mins,'moment_rows',momentRows,'seconds',toc(started),'scope','floating numerical diagnostics; interval and all-convex-body assertions are separately proved');
fid=fopen(fullfile(root,'results','quantitative_gap_matlab.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);disp(out);
