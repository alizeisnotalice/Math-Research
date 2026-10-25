% Fixed parity/signed-root interface check. No R4 union recomputation or search.
vg8_lo = [0,0,0,0;1,0,0,0;-1,0,0,0;-2,0,0,0;0,-1,0,0;...
          -2,-1,-1,-1;-3,0,0,0;-4,0,0,0;-2,-2,1,-1;2,0,0,0;...
          -1,0,0,0;0,0,0,0];
vg8_hi = [0,0,0,0;1,0,0,0;-1,0,0,0;-2,0,0,0;1,0,1,0;...
          -1,1,1,1;2,1,1,1;-2,0,0,0;3,-1,2,1;6,1,1,1;...
          1,2,2,2;2,2,2,2];
vg8_counts = zeros(12,2);
for vg8_i = 1:12
    vg8_E = floor(vg8_hi(vg8_i,1)/2)-floor((vg8_lo(vg8_i,1)-1)/2);
    vg8_O = vg8_hi(vg8_i,1)-vg8_lo(vg8_i,1)+1-vg8_E;
    vg8_T = prod(vg8_hi(vg8_i,2:4)-vg8_lo(vg8_i,2:4)+1);
    vg8_counts(vg8_i,:) = [vg8_E,vg8_O]*vg8_T;
    vg8_enum = [0,0];
    for vg8_k = vg8_lo(vg8_i,1):vg8_hi(vg8_i,1)
        vg8_q = floor(vg8_k/2);
        vg8_p = vg8_k-2*vg8_q;
        assert(vg8_p == 0 || vg8_p == 1);
        assert(2*vg8_q+vg8_p == vg8_k);
        vg8_enum(vg8_p+1) = vg8_enum(vg8_p+1)+vg8_T;
    end
    assert(all(vg8_enum == vg8_counts(vg8_i,:)));
end
assert(sum(vg8_counts(:)) == 391);
assert(all(vg8_counts(11,:) == [27,54]) && all(vg8_counts(12,:) == [54,27]));
vg8_eta = 1/64;
vg8_q = 4095/4096;
vg8_tau = 3*(1-vg8_q^4);
vg8_r0 = 3/2;
vg8_r1 = -1/8;
vg8_rtrue = (81-27*vg8_tau)^(1/4)/2;
vg8_remainder = vg8_rtrue-vg8_r0-vg8_r1*vg8_tau;
vg8_Br = 576/36481;
vg8_bound = vg8_Br*vg8_tau^2;
vg8_exact_remainder = -301940739/2251799813685248;
vg8_cross_error = abs(vg8_remainder-vg8_exact_remainder);
assert(0 < vg8_tau && vg8_tau < vg8_eta);
assert(vg8_remainder < 0 && abs(vg8_remainder) <= vg8_bound && vg8_cross_error < 1e-14);
assert((1+vg8_tau)+(1-vg8_tau) == 2);
disp('vg8 even/odd counts for the twelve fixed blocks');
disp(vg8_counts);
fprintf('vg8 atoms %g r1 %.17g tau %.17g remainder %.17g bound %.17g cross_error %.17g\n',...
    sum(vg8_counts(:)),vg8_r1,vg8_tau,vg8_remainder,vg8_bound,vg8_cross_error);
