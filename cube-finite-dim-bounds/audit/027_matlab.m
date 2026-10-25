% Fixed source-bound cross-check only; parent runs serially. No reset/path edits.
aud027T=4*log(2)-1;
aud027c=aud027T/(24*pi^3);
aud027R=15/(2*aud027T);
aud027q=sqrt(3*aud027R-2);
aud027a=(-1-aud027q)/3;aud027b=(-1+aud027q)/3;
aud027Im=integral(@(aud027u) 1-aud027R*aud027u.^2./(aud027u.^2-2*aud027u+3),0,-1/aud027a,'AbsTol',1e-13,'RelTol',1e-12);
aud027Ip=integral(@(aud027u) 1-aud027R*aud027u.^2./(aud027u.^2+2*aud027u+3),0,1/aud027b,'AbsTol',1e-13,'RelTol',1e-12);
aud027I=aud027Im+aud027Ip;
aud027Lb=4/pi*integral(@(aud027s) (2+3*aud027s).^2,-2/3,0,'AbsTol',1e-13,'RelTol',1e-12);
aud027mu=4*aud027c*aud027I/pi;
aud027cap=aud027mu/aud027Lb;
aud027negative=aud027mu-aud027Lb/200;
assert(abs(aud027Lb-32/(9*pi))<1e-12);
assert(aud027cap<.003347 && aud027negative<-.001871);
fprintf('AUD027 mu_upper=%.17g epsilon_cap=%.17g Lb=%.17g source_upper=%.17g PASS\n',aud027mu,aud027cap,aud027Lb,aud027negative);
