% Fixed checkerboard interface checks, not a new R4 certificate or a search.
vg11_L = [0,0;1,0;0,1;1,1;-1,-1;-2,-1;-1,-2;-2,-2;...
          -3,-2;-4,-3;-2,-1;2,-3;-1,-1;-1,0;0,-1;0,0;...
          -2,-2;-2,0;0,-2;-5,-4;-1,-1;-2,-2;-2,-1;1,0];
vg11_U = [0,0;1,0;0,1;1,1;-1,-1;-2,-1;-1,-2;-2,-2;...
          2,1;-2,-1;3,3;6,2;1,1;1,2;2,1;2,2;...
          -1,-1;-1,2;2,-1;-3,-2;0,0;2,2;2,1;0,0];
vg11_T = [ones(1,8),4,1,4,1,9,9,9,9,4,4,4,4,4,1,4,1];
vg11_cosets = [0,0;0,1;1,0;1,1];
vg11_signs = [1;-1;-1;1];
vg11_counts = zeros(24,4);
for vg11_i = 1:24
    for vg11_c = 1:4
        vg11_n = max(0,floor((vg11_U(vg11_i,:)-vg11_cosets(vg11_c,:))/2)-...
                      ceil((vg11_L(vg11_i,:)-vg11_cosets(vg11_c,:))/2)+1);
        vg11_counts(vg11_i,vg11_c) = prod(vg11_n)*vg11_T(vg11_i);
    end
    vg11_enum = zeros(1,4);
    for vg11_k1 = vg11_L(vg11_i,1):vg11_U(vg11_i,1)
        for vg11_k2 = vg11_L(vg11_i,2):vg11_U(vg11_i,2)
            vg11_p = mod(vg11_k1,2); vg11_q = mod(vg11_k2,2);
            assert(2*floor(vg11_k1/2)+vg11_p == vg11_k1);
            assert(2*floor(vg11_k2/2)+vg11_q == vg11_k2);
            vg11_enum(2*vg11_p+vg11_q+1) = vg11_enum(2*vg11_p+vg11_q+1)+vg11_T(vg11_i);
        end
    end
    assert(all(vg11_enum == vg11_counts(vg11_i,:)));
    for vg11_e = [-1/64,0,1/64]
        vg11_mass = vg11_counts(vg11_i,:)*(1+vg11_e*vg11_signs);
        assert(abs(vg11_mass-sum(vg11_counts(vg11_i,:))-vg11_e*vg11_counts(vg11_i,:)*vg11_signs) < 1e-12);
        assert(sum(1+vg11_e*vg11_signs) == 4);
    end
end
assert(sum(vg11_counts(:)) == 788);
assert(all(vg11_counts(14,:) == [18,9,36,18]));
vg11_root = 4095/4096;
vg11_tau = 9*(1-vg11_root^4);
vg11_r0 = 3/2; vg11_r1 = -1/24;
vg11_r = (81-9*vg11_tau)^(1/4)/2;
vg11_Br = 576/330625;
vg11_rem = vg11_r-vg11_r0-vg11_r1*vg11_tau;
assert(vg11_rem < 0 && abs(vg11_rem) <= vg11_Br*vg11_tau^2);
vg11_keylo = [0,0,0,0;0,0,0,0;-1,0,0,0;-1,0,0,0];
vg11_keyhi = [0,1,2,2;0,0,2,2;0,1,2,2;0,0,2,2];
vg11_Q = [2,2,1,1];
vg11_h = [1/8,-1/10,1/12,0];
vg11_u = vg11_h./vg11_Q;
vg11_phase0 = [vg11_cosets/2,zeros(4,2)];
vg11_phase = vg11_phase0+vg11_tau*vg11_signs*vg11_u;
vg11_lower = max(vg11_keyhi+vg11_phase,[],1)-vg11_r./vg11_Q;
vg11_upper = min(vg11_keylo+vg11_phase,[],1)+vg11_r./vg11_Q;
vg11_lower0 = max(vg11_keyhi+vg11_phase0,[],1)-vg11_r0./vg11_Q;
vg11_upper0 = min(vg11_keylo+vg11_phase0,[],1)+vg11_r0./vg11_Q;
vg11_lv = abs(vg11_u)-vg11_r1./vg11_Q;
vg11_uv = -abs(vg11_u)+vg11_r1./vg11_Q;
vg11_endpoint_error = max(abs([vg11_lower-vg11_lower0-vg11_tau*vg11_lv;...
                              vg11_upper-vg11_upper0-vg11_tau*vg11_uv]),[],1);
assert(all(vg11_endpoint_error <= vg11_Br./vg11_Q*vg11_tau^2+1e-14));
vg11_volume = prod(vg11_upper-vg11_lower);
vg11_expected = 3623344326421664984405288101222775868140421246247/...
               14615016373309029182036848327162830196559325429760;
vg11_cross_error = abs(vg11_volume-vg11_expected);
assert(vg11_cross_error < 1e-13);
vg11_anchored = vg11_phase-vg11_tau*ones(4,1)*vg11_u;
vg11_al = max(vg11_keyhi+vg11_anchored,[],1)-vg11_r./vg11_Q;
vg11_au = min(vg11_keylo+vg11_anchored,[],1)+vg11_r./vg11_Q;
assert(max(abs(vg11_al-vg11_lower+vg11_tau*vg11_u)) < 1e-14);
assert(max(abs(vg11_au-vg11_upper+vg11_tau*vg11_u)) < 1e-14);
disp('vg11 counts for the 24 fixed blocks');disp(vg11_counts);
fprintf('vg11 atoms %g r1 %.17g tau %.17g remainder %.17g bound %.17g\n',...
        sum(vg11_counts(:)),vg11_r1,vg11_tau,vg11_rem,vg11_Br*vg11_tau^2);
fprintf('vg11 receiver volume %.17g cross_error %.17g max_endpoint_error %.17g\n',...
        vg11_volume,vg11_cross_error,max(vg11_endpoint_error));
