% One rational simultaneous-event model; numerical crosscheck of exact identities.
% The finite quotients below check known affine formulas, not infer a derivative.
vg11b_a = 683/20000;
vg11b_b = 23/2000;
vg11b_base = [0,2;0,1];
vg11b_velocity = [vg11b_a,-vg11b_a;-vg11b_a-vg11b_b,-vg11b_a+vg11b_b];
vg11b_quotients = zeros(1,2);
vg11b_shielded = zeros(1,2);
vg11b_steps = [-1/100,1/100];
for vg11b_i = 1:2
    vg11b_dt = vg11b_steps(vg11b_i);
    vg11b_intervals = vg11b_base+vg11b_dt*vg11b_velocity;
    for vg11b_shield = 0:1
        vg11b_rows = vg11b_intervals;
        if vg11b_shield == 1
            vg11b_rows(3,:) = [-1,3];
        end
        vg11b_cuts = unique(vg11b_rows(:));
        vg11b_length = 0;
        for vg11b_j = 1:(numel(vg11b_cuts)-1)
            vg11b_mid = (vg11b_cuts(vg11b_j)+vg11b_cuts(vg11b_j+1))/2;
            if any(vg11b_rows(:,1) < vg11b_mid & vg11b_mid < vg11b_rows(:,2))
                vg11b_length = vg11b_length+vg11b_cuts(vg11b_j+1)-vg11b_cuts(vg11b_j);
            end
        end
        if vg11b_shield == 0
            vg11b_quotients(vg11b_i) = (vg11b_length-2)/vg11b_dt;
        else
            vg11b_shielded(vg11b_i) = (vg11b_length-4)/vg11b_dt;
        end
    end
end
vg11b_expected = [-2*vg11b_a,vg11b_b];
vg11b_error = max(abs(vg11b_quotients-vg11b_expected));
assert(vg11b_error < 1e-10 && max(abs(vg11b_shielded)) < 1e-10);
vg11b_F0upper = 1171/1000-1224/1000;
vg11b_F1lower = 1707/25000+1180/1000-1225/1000;
assert(vg11b_F0upper < 0 && vg11b_F1lower > 0);
fprintf('vg11b affine fixture left %.17g right %.17g max_cross_error %.17g shielded %.17g %.17g\n',...
    vg11b_quotients(1),vg11b_quotients(2),vg11b_error,vg11b_shielded(1),vg11b_shielded(2));
fprintf('vg11b actual-path cap signs F0upper %.17g F1lower %.17g; no threshold solve or path scan\n',...
    vg11b_F0upper,vg11b_F1lower);
