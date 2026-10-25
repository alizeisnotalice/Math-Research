% cb7_dual numerical all-cell verification only; no solver, clear, or path change.
cb7_dual_root=fileparts(mfilename('fullpath'));
cb7_dual_BN=[0 0 0 1;1 2 4 4;4 4 2 1;1 0 0 0];
cb7_dual_DN=[0 1 2 3;3 0 -3 -6;-6 -3 0 3;3 2 1 0];
for cb7_dual_q=1:2
    cb7_dual_rec=jsondecode(fileread(fullfile(cb7_dual_root,sprintf('cb7_dual_certificate_q%d.json',cb7_dual_q))));
    cb7_dual_K=cb7_dual_rec.K;
    cb7_dual_a=str2double(cb7_dual_rec.common_coefficient_numerators)/str2double(cb7_dual_rec.common_coefficient_denominator);
    cb7_dual_orbits=cb7_dual_rec.orbits;
    assert(all(cb7_dual_a>=0) && all(max(cb7_dual_orbits,[],2)>=2));
    cb7_dual_grid=zeros(2*cb7_dual_K+1,2*cb7_dual_K+1,2*cb7_dual_K+1);
    cb7_dual_origin=0;
    for cb7_dual_i=-cb7_dual_K:cb7_dual_K
        for cb7_dual_j=-cb7_dual_K:cb7_dual_K
            for cb7_dual_k=-cb7_dual_K:cb7_dual_K
                cb7_dual_v=[cb7_dual_i cb7_dual_j cb7_dual_k];
                if max(abs(cb7_dual_v))>=2
                    [cb7_dual_tf,cb7_dual_idx]=ismember(sort(abs(cb7_dual_v)),cb7_dual_orbits,'rows');
                    assert(cb7_dual_tf);
                    cb7_dual_grid(cb7_dual_i+cb7_dual_K+1,cb7_dual_j+cb7_dual_K+1,cb7_dual_k+cb7_dual_K+1)=cb7_dual_a(cb7_dual_idx);
                    cb7_dual_zero_factors=(2/3)*(cb7_dual_v==0)+(1/6)*(abs(cb7_dual_v)==1);
                    cb7_dual_origin=cb7_dual_origin+cb7_dual_a(cb7_dual_idx)*prod(cb7_dual_zero_factors);
                end
            end
        end
    end
    assert(cb7_dual_origin==0);
    cb7_dual_maximum=0;cb7_dual_integral=0;cb7_dual_cells=0;
    for cb7_dual_l1=-cb7_dual_K-2:cb7_dual_K+1
        for cb7_dual_l2=-cb7_dual_K-2:cb7_dual_K+1
            for cb7_dual_l3=-cb7_dual_K-2:cb7_dual_K+1
                cb7_dual_cell=[cb7_dual_l1 cb7_dual_l2 cb7_dual_l3];
                cb7_dual_poly=zeros(4,4,4);
                for cb7_dual_i=max(-cb7_dual_K,cb7_dual_l1-1):min(cb7_dual_K,cb7_dual_l1+2)
                    for cb7_dual_j=max(-cb7_dual_K,cb7_dual_l2-1):min(cb7_dual_K,cb7_dual_l2+2)
                        for cb7_dual_k=max(-cb7_dual_K,cb7_dual_l3-1):min(cb7_dual_K,cb7_dual_l3+2)
                            cb7_dual_coeff=cb7_dual_grid(cb7_dual_i+cb7_dual_K+1,cb7_dual_j+cb7_dual_K+1,cb7_dual_k+cb7_dual_K+1);
                            if cb7_dual_coeff==0, continue; end
                            cb7_dual_m=cb7_dual_cell-[cb7_dual_i cb7_dual_j cb7_dual_k]+3;
                            cb7_dual_bx=reshape(cb7_dual_BN(cb7_dual_m(1),:),[4 1 1]);
                            cb7_dual_by=reshape(cb7_dual_BN(cb7_dual_m(2),:),[1 4 1]);
                            cb7_dual_bz=reshape(cb7_dual_BN(cb7_dual_m(3),:),[1 1 4]);
                            cb7_dual_dx=reshape(cb7_dual_DN(cb7_dual_m(1),:),[4 1 1]);
                            cb7_dual_dy=reshape(cb7_dual_DN(cb7_dual_m(2),:),[1 4 1]);
                            cb7_dual_dz=reshape(cb7_dual_DN(cb7_dual_m(3),:),[1 1 4]);
                            cb7_dual_poly=cb7_dual_poly+(cb7_dual_coeff*cb7_dual_q^2/108)*(cb7_dual_dx.*cb7_dual_by.*cb7_dual_bz+cb7_dual_bx.*cb7_dual_dy.*cb7_dual_bz+cb7_dual_bx.*cb7_dual_by.*cb7_dual_dz);
                        end
                    end
                end
                cb7_dual_maximum=max(cb7_dual_maximum,max(cb7_dual_poly(:)));
                cb7_dual_cells=cb7_dual_cells+1;
                if all(cb7_dual_cell>=-cb7_dual_q & cb7_dual_cell<cb7_dual_q)
                    cb7_dual_integral=cb7_dual_integral+sum(cb7_dual_poly(:))/(64*cb7_dual_q^3);
                end
            end
        end
    end
    cb7_dual_value=1-cb7_dual_integral/8;
    cb7_dual_parts=split(string(cb7_dual_rec.D_exact),'/');
    cb7_dual_expected=str2double(cb7_dual_parts(1));
    if numel(cb7_dual_parts)==2, cb7_dual_expected=cb7_dual_expected/str2double(cb7_dual_parts(2)); end
    assert(cb7_dual_cells==(2*cb7_dual_K+4)^3 && cb7_dual_maximum<=1+1e-10);
    assert(abs(cb7_dual_value-cb7_dual_expected)<1e-10);
    fprintf('cb7_dual q=%d PASS: %d full cells, max Bernstein %.17g, psi0 %.17g, D %.17g\n',cb7_dual_q,cb7_dual_cells,cb7_dual_maximum,cb7_dual_origin,cb7_dual_value);
end
% Floating cross-check only. Exact validity is from the rational certificate/checker.
