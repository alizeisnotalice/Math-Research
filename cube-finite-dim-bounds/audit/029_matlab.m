% Parent serial finite arithmetic crosscheck. No candidate/source cover or LP.
% Can use eval(fileread(...)); no local functions, clear, or path mutation.
aud029_V=sym([1 -1 -1;1 -1 1;1 1 -1;1 1 1]).';
aud029_triples=nchoosek(1:4,3);
for aud029_i=1:4
    aud029_M=aud029_V(:,aud029_triples(aud029_i,:));
    aud029_inverse=inv(aud029_M);
    assert(abs(det(aud029_M))==4);
    assert(all(isAlways(sum(aud029_inverse.^2,2)<=1)));
    aud029_remaining=setdiff(1:4,aud029_triples(aud029_i,:));
    assert(all(isAlways(abs(aud029_inverse*aud029_V(:,aud029_remaining))==1)));
end
aud029_r=sym('aud029_r','real');
aud029_m0=int(4*sym(pi)*aud029_r^2*(1-aud029_r^2)^2,aud029_r,0,1);
aud029_m1=6*sym(1)/8*aud029_m0;
assert(isAlways(aud029_m0==32*sym(pi)/105));
assert(isAlways(aud029_m1==8*sym(pi)/35));
aud029_tau=sym('aud029_tau','real');
aud029_p=sym([1;1;1])+aud029_tau*sym([1;1;-1]);
aud029_dist2=sym(zeros(1,4));
for aud029_i=1:4
    aud029_dist2(aud029_i)=expand(aud029_p.'*aud029_p- ...
        (aud029_p.'*aud029_V(:,aud029_i))^2/3);
end
assert(isAlways(aud029_dist2(1)==8*(aud029_tau^2+aud029_tau+1)/3));
assert(isAlways(aud029_dist2(2)==aud029_dist2(1)));
assert(isAlways(aud029_dist2(3)==sym(8)/3));
assert(isAlways(aud029_dist2(4)==8*aud029_tau^2/3));
aud029_far_coarse= (4*sym(2)/3-1)/48 -48*(sym(5)/2)/60^2;
assert(isAlways(aud029_far_coarse==sym(1)/720));
aud029_c=(4*log(2)-1)/(24*pi^3);
aud029_tailgap=aud029_c/2-48*sqrt(6)/(pi^3*60^2);
assert(aud029_tailgap>1/(720*pi^3));
aud029_supportgap=5/(16*pi^3)/(25/16)-aud029_c;
aud029_buffergap=5/(16*pi^3)/(9/4)-aud029_c;
assert(aud029_supportgap>3/(40*pi^3));
assert(aud029_buffergap>1/(72*pi^3));
aud029_LBcoefficient=4*(2*(sym(5)/4)/sqrt(sym(3)))*4/sym(pi);
assert(isAlways(aud029_LBcoefficient==40/(sqrt(sym(3))*sym(pi))));
fprintf('029 PASS all4bases, masses32pi/105 and8pi/35, transverse min8/3, coarse gap1/720\n');
fprintf('supportgap=%.17g buffergap=%.17g tailgap_t60=%.17g\n', ...
    aud029_supportgap,aud029_buffergap,aud029_tailgap);
