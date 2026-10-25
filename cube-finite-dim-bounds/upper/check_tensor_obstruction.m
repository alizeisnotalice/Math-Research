% U01 independent floating-point check; all variables use shared-session up prefix.
upN = [3,4];
upMass = -(5/16).^(upN-1);
upS = pi/2;
upT = pi;
upQs = (sin(upS)/upS)^2;
upQt = (sin(upT)/upT)^2;
upFs0 = upS^2*upQs;
upF0t = upT^2*upQt;
upFst = (upS^2+upT^2)*upQs*upQt;
fprintf('U01 n=3,4 tent mass: %.0f %.0f\n',2.^upN);
fprintf('U01 negative slab masses: %.14f %.14f\n',upMass);
fprintf('U01 Fourier Fs0=%.16g F0t=%.16g Fst=%.16g defect=%.16g\n', ...
    upFs0,upF0t,upFst,abs(upFst-upFs0)-upF0t);
assert(all(upMass<0));
assert(abs(upFst-upFs0)-upF0t>0.999999);
