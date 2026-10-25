% Serialized by primary agent; no clear/path/cwd modification.
n4r6a = [1+sqrt(2)/2,(2+2^(3/4)*3^(1/4))/4,(1+2^(1/4))/2,(1+2^(1/4))/2];
n4r6source = [1,4,16,24];
n4r6raw = mod(n4r6source.^(1/4)./n4r6a,1);
n4r6phi = min(n4r6raw,1-n4r6raw);
n4r6eps = 1e-4; n4r6singletonr = n4r6eps^(1/4)/2;
n4r6singletonVolume = (2*n4r6singletonr)^4;
assert(abs(n4r6singletonVolume-n4r6eps)<1e-15);
assert(all(n4r6phi>=0 & n4r6phi<=.5));
fprintf('n4 round006 mixed phi=[%.16g %.16g %.16g %.16g], singleton volume=%.16g, mass denominator=%.16g\n',n4r6phi,n4r6singletonVolume,1+n4r6eps);
