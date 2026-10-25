% cb7 finite floating cross-check only; no clear/path changes, all cb7 names.
cb7_nlist=[3,4];
cb7_A=[0.4089479272201547,0.18157763153667977];
cb7_B4=[0.005004477192161995,0.004];
cb7_B6=[0.0009386738283819751,0.0005];
cb7_R=[16/5,18/5];
for cb7_idx=1:2
    cb7_n=cb7_nlist(cb7_idx);cb7_rad=cb7_R(cb7_idx);
    cb7_d=cb7_B6(cb7_idx)*cb7_rad^6;
    cb7_e=cb7_B4(cb7_idx)*cb7_rad^4-15*cb7_d/(cb7_n+8);
    cb7_C=-3*cb7_B4(cb7_idx)*cb7_rad^4/(cb7_n+2)+30*cb7_d/((cb7_n+8)*(cb7_n+4));
    cb7_values=[];
    for cb7_s=1:cb7_n
        cb7_values(end+1)=cb7_d/cb7_s^2+cb7_e/cb7_s+cb7_C;
        cb7_T=-2*cb7_e/(3*cb7_d);
        for cb7_p=1:cb7_s-1
            cb7_q=cb7_s-cb7_p;
            if cb7_p~=cb7_q
                cb7_u=(1-cb7_q*cb7_T)/(cb7_p-cb7_q);cb7_v=cb7_T-cb7_u;
                if cb7_u>0 && cb7_v>0
                    cb7_values(end+1)=cb7_d*(cb7_p*cb7_u^3+cb7_q*cb7_v^3)+cb7_e*(cb7_p*cb7_u^2+cb7_q*cb7_v^2)+cb7_C;
                end
            end
        end
    end
    cb7_margin=cb7_A(cb7_idx)-cb7_rad^(2-cb7_n)+min(cb7_values);
    assert(cb7_margin>0);
    fprintf('cb7 n=%d outer-sphere margin %.17g (%d candidate entries)\n',cb7_n,cb7_margin,numel(cb7_values));
end
cb7_poly3=[-1/196,0,23/420,0,0,0,-1543/1470,1];
cb7_poly4=[-1/11,25/66,0,-85/66,1];
cb7_roots3=roots(cb7_poly3);cb7_roots4=roots(cb7_poly4);
cb7_pos3=sort(real(cb7_roots3(abs(imag(cb7_roots3))<1e-8 & real(cb7_roots3)>0)));
cb7_pos4=sort(real(cb7_roots4(abs(imag(cb7_roots4))<1e-8 & real(cb7_roots4)>0)));
assert(numel(cb7_pos3)==3 && max(abs(cb7_pos3-[1;2;3]))<1e-8);
assert(numel(cb7_pos4)==3 && max(abs(cb7_pos4-[1;2;3]))<1e-8);
cb7_t=[1/4,3/4];cb7_cube_reduced=(4*cb7_t-3).*(4*cb7_t-1)/12;
assert(all(cb7_cube_reduced==0));
fprintf('cb7 PASS: no-guard radial examples each have positive roots [1,2,3]; cube equal-multiplicity roots [1/4,3/4].\n');
% These checks do not certify cube global minima or angular mass/cost.
