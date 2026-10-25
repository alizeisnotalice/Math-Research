% Fixed dual-satellite local toy. Numerical crosscheck, not a new lower bound.
vg7_weights = [1/2,1/2];
vg7_rprime = sum(vg7_weights)/8;
vg7_d = 3/16;
vg7_Bsingle = 1/16-vg7_d/2;
vg7_Bjoint = vg7_rprime-vg7_d/2;
vg7_L = -3/32; vg7_R = 3/32;
vg7_D0 = -(vg7_R-vg7_L);
vg7_D1 = vg7_D0+2*max(0,vg7_Bsingle);
vg7_D12 = vg7_D0+2*max(0,vg7_Bjoint);
vg7_mixed = vg7_D12-2*vg7_D1+vg7_D0;
vg7_area = 1/8; vg7_V0 = 1/2;
vg7_G0 = vg7_area*vg7_D0+vg7_V0*vg7_d;
vg7_G1 = vg7_area*vg7_D1+vg7_V0*(vg7_d-1/2);
vg7_G12 = vg7_area*vg7_D12+vg7_V0*(vg7_d-1);
vg7_q = 2049/2048;
vg7_sigma = vg7_q^4-1;
vg7_h = 1+vg7_d*vg7_sigma;
vg7_front = 1/(2*vg7_h)-1/2;
vg7_born = vg7_q/(2*vg7_h)-1/2;
vg7_intervals = [-2,vg7_front;-vg7_front,2;-vg7_born,vg7_born];
vg7_cuts = unique(vg7_intervals(:));
vg7_length = 0;
for vg7_i = 1:(numel(vg7_cuts)-1)
    vg7_mid = (vg7_cuts(vg7_i)+vg7_cuts(vg7_i+1))/2;
    if any(vg7_intervals(:,1) < vg7_mid & vg7_mid < vg7_intervals(:,2))
        vg7_length = vg7_length+vg7_cuts(vg7_i+1)-vg7_cuts(vg7_i);
    end
end
vg7_volume = vg7_area*vg7_length;
vg7_expected = 563121827635203/1126312525791244;
vg7_cross_error = abs(vg7_volume-vg7_expected);
assert(vg7_Bsingle < 0 && vg7_Bjoint == 1/32);
assert(vg7_mixed == 1/16 && vg7_G12-2*vg7_G1+vg7_G0 == 1/128);
assert(0 < vg7_born && vg7_born < -vg7_front && vg7_cross_error < 1e-12);
fprintf('vg7 rprime %.17g Bsingle %.17g Bjoint %.17g mixed %.17g quotient_mixed %.17g\n',...
         vg7_rprime,vg7_Bsingle,vg7_Bjoint,vg7_mixed,vg7_G12-2*vg7_G1+vg7_G0);
fprintf('vg7 sigma %.17g volume %.17g cross_error %.17g\n',...
         vg7_sigma,vg7_volume,vg7_cross_error);
