% U06 independent all-root radial integration; main agent runs MATLAB serially.
upU06A=.4089479272201547;
upU06B4=.005004477192161995;
upU06B6=.0009386738283819751;
upU06R=16/5;
upU06mass=8*integral2(@(upU06u,upU06phi) upU06radial(upU06u,upU06phi,upU06A,upU06B4,upU06B6,upU06R), ...
    0,1,0,pi/2,'AbsTol',1e-9,'RelTol',3e-8);
upU06H4corner=3-3*9/5;
upU06H6corner=3-(15/11)*9+(30/77)*27;
upU06corner=1/sqrt(3)-upU06A-upU06B4*upU06H4corner-upU06B6*upU06H6corner;
fprintf('U06 n3 mass_dimless=%.17g corner_value_dimless=%.17g corner_ratio=%.17g\n',upU06mass,upU06corner,upU06mass/(8*upU06corner));
[upU06checks,upU06counts]=upU06radial([1,0,1/sqrt(3),.37],[0,pi/4,pi/4,.29],upU06A,upU06B4,upU06B6,upU06R);
fprintf('U06 axis/edge/diagonal/generic radial positive interval counts: %d %d %d %d\n',upU06counts);

function [upU06v,upU06counts]=upU06radial(upU06u,upU06phi,upU06A,upU06B4,upU06B6,upU06R)
upU06y1=(1-upU06u.^2).*cos(upU06phi).^2;
upU06y2=(1-upU06u.^2).*sin(upU06phi).^2;
upU06y3=upU06u.^2;
upU06s2=upU06y1.^2+upU06y2.^2+upU06y3.^2;
upU06h4=upU06s2-3/5;
upU06h6=upU06y1.^3+upU06y2.^3+upU06y3.^3-(15/11)*upU06s2+30/77;
upU06v=zeros(size(upU06h4));upU06counts=zeros(size(upU06h4));
for upU06j=1:numel(upU06h4)
    upU06coef=[-upU06B6*upU06h6(upU06j)*upU06R^7,0,-upU06B4*upU06h4(upU06j)*upU06R^5,0,0,0,-upU06A*upU06R,1];
    upU06coef=upU06coef(find(upU06coef~=0,1):end);
    upU06roots=roots(upU06coef);
    upU06roots=sort(real(upU06roots(abs(imag(upU06roots))<1e-7 & real(upU06roots)>0 & real(upU06roots)<1)));
    upU06cuts=[0;upU06roots;1];
    for upU06k=1:numel(upU06cuts)-1
        if polyval(upU06coef,(upU06cuts(upU06k)+upU06cuts(upU06k+1))/2)>0
            upU06lo=upU06R*upU06cuts(upU06k);upU06hi=upU06R*upU06cuts(upU06k+1);
            upU06J=@(upU06r) upU06r.^2/2-upU06A*upU06r.^3/3-upU06B4*upU06h4(upU06j)*upU06r.^7/7-upU06B6*upU06h6(upU06j)*upU06r.^9/9;
            upU06v(upU06j)=upU06v(upU06j)+upU06J(upU06hi)-upU06J(upU06lo);
            upU06counts(upU06j)=upU06counts(upU06j)+1;
        end
    end
end
end
