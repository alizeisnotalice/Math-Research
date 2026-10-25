% Parent serial real-Hermite DD identity/limit cross only, not a spatial grid.
% No local functions, clear/path changes, Green module or cost evaluation.
aud032_e=sym('aud032_e','positive');
aud032_b=sym('aud032_b',[1,3],'real');
aud032_derivative=sym(0);
for aud032_j=1:3
    aud032_other=setdiff(1:3,aud032_j);
    aud032_derivative=aud032_derivative+aud032_e/(aud032_e^2+aud032_b(aud032_j)^2)/ ...
        prod(aud032_b(aud032_j)-aud032_b(aud032_other));
end
aud032_pairs=aud032_b(1)*aud032_b(2)+aud032_b(1)*aud032_b(3)+aud032_b(2)*aud032_b(3);
aud032_realproduct=aud032_e*(aud032_e^2-aud032_pairs)/prod(aud032_e^2+aud032_b.^2);
assert(isAlways(simplify(aud032_derivative+aud032_realproduct)==0));
aud032_examples=sym([0 0 2;0 0 -2;0 2 2;0 -2 -2;2 2 2;-2 -2 -2;0 2 -3;2 2 -3]);
aud032_gamma=[-sym(1)/4;-sym(1)/4;sym(1)/4;sym(1)/4;0;0;-sym(1)/6;0];
aud032_finite=[log(sym(2))/4;log(sym(2))/4;(1-log(sym(2)))/4; ...
    (1-log(sym(2)))/4;-sym(1)/8;-sym(1)/8; ...
    log(sym(2))/10+log(sym(3))/15;(log(sym(3))-log(sym(2))+sym(5)/2)/25];
for aud032_case=1:size(aud032_examples,1)
    aud032_nodes=sort(aud032_examples(aud032_case,:));aud032_DD=sym(zeros(3));
    for aud032_i=1:3
        aud032_DD(aud032_i,1)=log(aud032_e^2+aud032_nodes(aud032_i)^2)/2;
    end
    for aud032_order=1:2
        for aud032_i=1:3-aud032_order
            aud032_j=aud032_i+aud032_order;
            if aud032_nodes(aud032_i)==aud032_nodes(aud032_j)
                aud032_node=aud032_nodes(aud032_i);
                if aud032_order==1
                    aud032_DD(aud032_i,2)=aud032_node/(aud032_e^2+aud032_node^2);
                else
                    aud032_DD(aud032_i,3)=(aud032_e^2-aud032_node^2)/ ...
                        (2*(aud032_e^2+aud032_node^2)^2);
                end
            else
                aud032_DD(aud032_i,aud032_order+1)= ...
                    (aud032_DD(aud032_i+1,aud032_order)-aud032_DD(aud032_i,aud032_order))/ ...
                    (aud032_nodes(aud032_j)-aud032_nodes(aud032_i));
            end
        end
    end
    aud032_limit=limit(aud032_DD(1,3)-aud032_gamma(aud032_case)*log(aud032_e),aud032_e,0,'right');
    assert(isAlways(simplify(aud032_limit-aud032_finite(aud032_case))==0));
end
aud032_tau=sym(1)/3;
aud032_peak=5/(8*(1+aud032_tau)*(2+3*aud032_tau)*(4+5*aud032_tau));
assert(isAlways(aud032_peak==sym(15)/544));
fprintf('032 PASS real rational DD derivative;8 signed scalar confluence limits;peak_tau1_3=15/544\n');
