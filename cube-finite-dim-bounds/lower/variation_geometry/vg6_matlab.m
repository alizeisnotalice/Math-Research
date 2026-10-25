% One fixed toy finite-epsilon certificate. Numerical crosscheck only.
vg6_e = 1/256;
vg6_s = 1/4;
vg6_r = vg6_s/2;
vg6_lo = [-2,0,0,0;2*vg6_e+vg6_e^2,0,0,0;-vg6_e/2,0,0,0];
vg6_hi = [-vg6_e,1+vg6_e+vg6_e^2,1+2*vg6_e,1+3*vg6_e;...
          2,1-vg6_e,1+vg6_e,1-2*vg6_e;vg6_e/2,1,1,1];
vg6_center = [-2,0,0,0]+vg6_e*[1/3,-1/4,1/5,-1/6];
vg6_lo(4,:) = vg6_center-vg6_r;
vg6_hi(4,:) = vg6_center+vg6_r;
vg6_volume = 0;
for vg6_mask = 1:15
    vg6_ids = find(bitget(vg6_mask,1:4));
    vg6_widths = max(0,min(vg6_hi(vg6_ids,:),[],1)-max(vg6_lo(vg6_ids,:),[],1));
    vg6_volume = vg6_volume+(-1)^(numel(vg6_ids)+1)*prod(vg6_widths);
end
vg6_C = 557311/2048;
vg6_K = 641/40;
vg6_D = 111/16;
vg6_error = abs(vg6_volume-4-vg6_D*vg6_e);
vg6_bound = vg6_C*vg6_e^2+vg6_K*vg6_s^7;
vg6_gain = vg6_volume/(1+vg6_e)-4;
vg6_guaranteed = 4257541/673710080;
vg6_expected = 24908410373519/6184752906240;
vg6_cross_error = abs(vg6_volume-vg6_expected);
assert(vg6_error <= vg6_bound && vg6_gain > vg6_guaranteed && vg6_guaranteed > 0);
assert(vg6_cross_error < 1e-12);
fprintf('vg6 epsilon %.17g volume %.17g gain %.17g guaranteed %.17g\n',...
         vg6_e,vg6_volume,vg6_gain,vg6_guaranteed);
fprintf('vg6 error %.17g bound %.17g cross_error %.17g\n',...
         vg6_error,vg6_bound,vg6_cross_error);
