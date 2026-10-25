% MATLAB numerical crosscheck only. No clear, path mutation, or symbolic toolbox.
vg5_t = [1/8,1/16,1/32,1/64];
vg5_a = zeros(size(vg5_t));
vg5_b = zeros(size(vg5_t));
vg5_exact_a = [31/4,63/8,127/16,255/32];
vg5_exact_b = [10327/2048,145715/32768,2204779/524288,34377947/8388608];
for vg5_i = 1:numel(vg5_t)
    vg5_e = vg5_t(vg5_i);
    vg5_intervals = [-4,-vg5_e; 2*vg5_e,4; -vg5_e/2,0; vg5_e,3*vg5_e/2];
    vg5_cuts = unique(vg5_intervals(:));
    for vg5_j = 1:(numel(vg5_cuts)-1)
        vg5_mid = (vg5_cuts(vg5_j)+vg5_cuts(vg5_j+1))/2;
        if any(vg5_intervals(:,1) < vg5_mid & vg5_mid < vg5_intervals(:,2))
            vg5_a(vg5_i) = vg5_a(vg5_i)+vg5_cuts(vg5_j+1)-vg5_cuts(vg5_j);
        end
    end
    vg5_lo = [-2,0,0,0;2*vg5_e,0,0,0;-vg5_e/2,0,0,0];
    vg5_hi = [-vg5_e,1+vg5_e,1+2*vg5_e,1+3*vg5_e;...
              2,1-vg5_e,1+vg5_e,1-2*vg5_e;vg5_e/2,1,1,1];
    for vg5_mask = 1:7
        vg5_ids = find(bitget(vg5_mask,1:3));
        vg5_widths = max(0,min(vg5_hi(vg5_ids,:),[],1)-max(vg5_lo(vg5_ids,:),[],1));
        vg5_b(vg5_i) = vg5_b(vg5_i)+(-1)^(numel(vg5_ids)+1)*prod(vg5_widths);
    end
end
vg5_error_a = max(abs(vg5_a-vg5_exact_a));
vg5_error_b = max(abs(vg5_b-vg5_exact_b));
assert(vg5_error_a < 1e-12 && vg5_error_b < 1e-12);
disp('vg5 epsilon, group A union length, group B 4D union volume');
disp([vg5_t(:),vg5_a(:),vg5_b(:)]);
fprintf('vg5 max errors A %.17g B %.17g\n',vg5_error_a,vg5_error_b);
