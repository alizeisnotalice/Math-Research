% U07 selected n3/n4, main-serialized MATLAB; guards/cube monotonicity exact in checker.
upU07params=[3,.41061707799980696,.0052535036230922875,.0009893885508363638,3.1643099997673794; ...
    4,.14344247730484339,.002778738033100787,.0003790027706453,3.3162574578337103];
for upU07i=1:2
    upU07n=upU07params(upU07i,1);upU07A=upU07params(upU07i,2);
    upU07B4=upU07params(upU07i,3);upU07B6=upU07params(upU07i,4);upU07R=upU07params(upU07i,5);
    if upU07n==3
        upU07mass=8*integral2(@(upU07theta,upU07phi) upU07radial(3,upU07theta,zeros(size(upU07theta)),upU07phi,upU07A,upU07B4,upU07B6,upU07R).*sin(upU07theta), ...
            0,pi/2,0,pi/2,'AbsTol',1e-9,'RelTol',1e-8);
    else
        upU07mass=16*integral3(@(upU07theta,upU07eta,upU07phi) upU07radial(4,upU07theta,upU07eta,upU07phi,upU07A,upU07B4,upU07B6,upU07R).*sin(upU07theta).^2.*sin(upU07eta), ...
            0,pi/2,0,pi/2,0,pi/2,'AbsTol',1e-8,'RelTol',1e-7);
    end
    upU07h4corner=upU07n-3*upU07n^2/(upU07n+2);
    upU07h6corner=upU07n-15*upU07n^2/(upU07n+8)+30*upU07n^3/((upU07n+8)*(upU07n+4));
    upU07cube=upU07n^((2-upU07n)/2)-upU07A-upU07B4*upU07h4corner-upU07B6*upU07h6corner;
    fprintf('U07 n=%d R=%.17g mass_dimless=%.17g cube_min_dimless=%.17g cost=%.17g\n',upU07n,upU07R,upU07mass,upU07cube,upU07mass/(2^upU07n*upU07cube));
end

function upU07v=upU07radial(upU07n,upU07theta,upU07eta,upU07phi,upU07A,upU07B4,upU07B6,upU07R)
if upU07n==3
    upU07y={sin(upU07theta).^2.*cos(upU07phi).^2,sin(upU07theta).^2.*sin(upU07phi).^2,cos(upU07theta).^2};
else
    upU07y={sin(upU07theta).^2.*sin(upU07eta).^2.*cos(upU07phi).^2,sin(upU07theta).^2.*sin(upU07eta).^2.*sin(upU07phi).^2,sin(upU07theta).^2.*cos(upU07eta).^2,cos(upU07theta).^2};
end
upU07s2=zeros(size(upU07theta));upU07s3=upU07s2;
for upU07j=1:upU07n
    upU07s2=upU07s2+upU07y{upU07j}.^2;upU07s3=upU07s3+upU07y{upU07j}.^3;
end
upU07h4=upU07s2-3/(upU07n+2);
upU07h6=upU07s3-15*upU07s2/(upU07n+8)+30/((upU07n+8)*(upU07n+4));
upU07lo=zeros(size(upU07h4));upU07hi=upU07R*ones(size(upU07h4));
for upU07j=1:60
    upU07mid=(upU07lo+upU07hi)/2;
    upU07positive=upU07mid.^(2-upU07n)-upU07A-upU07B4*upU07h4.*upU07mid.^4-upU07B6*upU07h6.*upU07mid.^6>0;
    upU07lo(upU07positive)=upU07mid(upU07positive);upU07hi(~upU07positive)=upU07mid(~upU07positive);
end
upU07r=(upU07lo+upU07hi)/2;
upU07v=upU07r.^2/2-upU07A*upU07r.^upU07n/upU07n-upU07B4*upU07h4.*upU07r.^(upU07n+4)/(upU07n+4)-upU07B6*upU07h6.*upU07r.^(upU07n+6)/(upU07n+6);
end
