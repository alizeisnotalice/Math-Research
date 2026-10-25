% 037 narrow exact scalar crosscheck; no clear/path or new integral/Green.
aud037I=sym(92072351)/sym(500000000000);
aud037U=sym(14371)/1000;
aud037p=sym(355)/113;
aud037B=256*aud037I^3*(sym(640)/17*aud037p^3)^4;
aud037gap=aud037U^3-aud037B;
assert(isAlways(aud037gap>0));
aud037ratios=[sym(1)/24,sym(17)/640,sym(482191)/13088544];
assert(isAlways(min(aud037ratios)==sym(17)/640));
aud037angle=-sym(1)/32*(-sym(17)/20);
assert(isAlways(aud037angle==sym(17)/640));
aud037a5=sym(0);aud037a239=sym(0);
for aud037k=0:12
 aud037a5=aud037a5+(-1)^aud037k/sym(5)^(2*aud037k+1)/(2*aud037k+1);
end
for aud037k=0:3
 aud037a239=aud037a239+(-1)^aud037k/sym(239)^(2*aud037k+1)/(2*aud037k+1);
end
assert(isAlways(16*aud037a5-4*aud037a239<aud037p));
fprintf('037 PASS minimum ratio17/640, anchor scalar, Machin; cube upper=%s gap=%s\n',char(aud037B),char(aud037gap));
disp(vpa((256*aud037I^3*(sym(640)/17*sym(pi)^3)^4)^(sym(1)/3),50));
