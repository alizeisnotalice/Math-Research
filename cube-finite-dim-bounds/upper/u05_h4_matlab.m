% U05 independent adaptive sphere integral. Main agent runs MATLAB serially.
% Actual F=(4*pi)^(-1)*(r^(-1)-A-BH4), W=1; common factor cancels in cost.
upU05A=0.39660079267158965;
upU05B=.003;
upU05R=(1/(4*(4/15)*upU05B))^(1/5);
upU05mass=8*integral2(@(upU05u,upU05phi) upU05radial(upU05u,upU05phi,upU05A,upU05B,upU05R), ...
    0,1,0,pi/2,'AbsTol',1e-10,'RelTol',1e-9);
upU05min=1/sqrt(3)+(12/5)*upU05B-upU05A;
upU05cost=upU05mass/(8*upU05min);
upU05outermax=1/upU05R-upU05A+(4/15)*upU05B*upU05R^4;
fprintf('U05 n3 A=%.17g B=%.17g R=%.17g mass_dimless=%.17g cube_min_dimless=%.17g cost=%.17g outer_max_dimless=%.17g\n', ...
    upU05A,upU05B,upU05R,upU05mass,upU05min,upU05cost,upU05outermax);

function upU05v=upU05radial(upU05u,upU05phi,upU05A,upU05B,upU05R)
upU05h=(1-upU05u.^2).^2.*(cos(upU05phi).^4+sin(upU05phi).^4)+upU05u.^4-3/5;
upU05lo=zeros(size(upU05h));
upU05hi=upU05R*ones(size(upU05h));
for upU05iteration=1:60
    upU05mid=(upU05lo+upU05hi)/2;
    upU05positive=1./upU05mid-upU05A-upU05B*upU05h.*upU05mid.^4>0;
    upU05lo(upU05positive)=upU05mid(upU05positive);
    upU05hi(~upU05positive)=upU05mid(~upU05positive);
end
upU05rho=(upU05lo+upU05hi)/2;
upU05v=upU05rho.^2/2-upU05A*upU05rho.^3/3-upU05B*upU05h.*upU05rho.^7/7;
end
