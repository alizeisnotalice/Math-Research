% Parent-serialized scalar cross-check of the independent moving-cube formulas.
aud10T=[1/8,1/16,1/32,1/64];aud10Rows=zeros(4,5);
for aud10I=1:4
 aud10R=aud10T(aud10I);aud10E=(2*aud10R)^4;
 aud10P=aud10E*[1/3,-1/4,1/5,-1/6];
 aud10Net=aud10E-prod(aud10R+aud10P-aud10E);
 aud10Error=aud10Net-15*aud10E/16;
 aud10Bound=233*aud10E*aud10R^3/60;
 aud10Overlap=16*aud10E*aud10R^3;
 assert(aud10Error>=0 && aud10Error<=aud10Bound*(1+1e-10));
 aud10Rows(aud10I,:)=[aud10E,aud10Net,aud10Error,aud10Bound,aud10Overlap];
end
fprintf('AUDIT010 columns epsilon, uncovered, error, bound, thin overlap\n');disp(aud10Rows);
aud10E=(2/16)^4;aud10Clipped=15*aud10E/16;aud10Periodic=0;aud10J=210;aud10Normalized=aud10Clipped/aud10J;
assert(abs(aud10J*aud10Normalized-aud10Clipped)<1e-18);
fprintf('AUDIT010 clipped=%0.17g full_periodic=%g normalized=%0.17g physical=%0.17g PASS\n',aud10Clipped,aud10Periodic,aud10Normalized,aud10J*aud10Normalized);
