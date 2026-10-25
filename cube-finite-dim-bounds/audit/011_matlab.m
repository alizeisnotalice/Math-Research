% Parent-serialized numerical cross-check of exact radial remainder identities.
aud11T=[1/4,1/2,3/4,9/10];aud11D3=3*(1-aud11T.^2);
aud11Left3=1./aud11T-1-(4*aud11D3+aud11D3.^2)/24;
aud11Right3=(1-aud11T).^3.*(3*aud11T.^2+9*aud11T+8)./(8*aud11T);
aud11D4=[1/4,1,2,7/2];
aud11Left4=1./(4-aud11D4)-1/4-(4*aud11D4+aud11D4.^2)/64;
aud11Right4=aud11D4.^3./(64*(4-aud11D4));
aud11Error3=max(abs(aud11Left3-aud11Right3));aud11Error4=max(abs(aud11Left4-aud11Right4));
assert(aud11Error3<1e-12 && aud11Error4<1e-12);
assert(all(aud11Right3>=0) && all(aud11Right4>=0));
fprintf('AUDIT011 radial-remainder identities PASS n3_error=%.17g n4_error=%.17g Bcrit3=%.17g Bcrit4=%.17g\n',aud11Error3,aud11Error4,5*sqrt(3)/144,1/32);
