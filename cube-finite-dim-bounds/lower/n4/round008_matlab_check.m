% Serialized by primary; n4r8-only variables, no clear/path/cwd.
n4r8P = [1+sqrt(2)/2,(2+2^(3/4)*3^(1/4))/4,(1+2^(1/4))/2,(1+2^(1/4))/2];
n4r8Q = n4r8P.*[2 1 1 1]; n4r8h=[.02 .02 .02 0]; n4r8d=n4r8P.*n4r8h;
n4r8phase=[.5+n4r8h(1)/2,n4r8h(2:4)];
assert(max(abs(n4r8Q.*n4r8phase-([n4r8P(1) 0 0 0]+n4r8d)))<1e-12);
n4r8L=-3; n4r8U=1; n4r8even=floor(n4r8U/2)-floor((n4r8L-1)/2); n4r8odd=n4r8U-n4r8L+1-n4r8even;
n4r8counts=6*[n4r8even n4r8odd]; n4r8eps=.05;
n4r8mass=n4r8counts*[1+n4r8eps;1-n4r8eps]; n4r8r=n4r8mass^(1/4)/2;
assert(n4r8even==2 && n4r8odd==3 && abs(n4r8mass-29.7)<1e-12);
fprintf('n4 round008 negative/cross0 parity=[%d %d], totalmass=%.16g r=%.16g supercellmass=%.16g physicald=[%.16g %.16g %.16g %.16g]\n',n4r8counts,n4r8mass,n4r8r,(1+n4r8eps)+(1-n4r8eps),n4r8d);
