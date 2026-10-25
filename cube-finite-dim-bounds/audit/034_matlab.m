% Parent serial theory-constant/chord/mean cross only; no live I parts or nodes.
% No local functions, clear/path mutation, cache or Green evaluation.
for aud034_n=[3,4]
    aud034_q=sym(aud034_n)/(aud034_n-1);
    aud034_exponent=((2-aud034_n)*(aud034_q-1)+1)/aud034_q;
    assert(isAlways(aud034_exponent==sym(1)/aud034_n));
    aud034_orthant=-(aud034_n-1)/aud034_q;
    assert(isAlways(aud034_orthant==-sym((aud034_n-1)^2)/aud034_n));
end
aud034_a=sym('aud034_a','positive');aud034_b=sym('aud034_b','positive');
aud034_identity=(aud034_a^2+aud034_a*aud034_b+aud034_b^2)/(aud034_a+aud034_b) ...
    -(aud034_a^3-aud034_b^3)/(aud034_a^2-aud034_b^2);
assert(isAlways(simplify(aud034_identity)==0));
aud034_u=sym('aud034_u','real');
for aud034_k=[0,1,2,32]
    aud034_mean=int(chebyshevT(aud034_k,2*aud034_u-1),aud034_u,0,1);
    if aud034_k==0,aud034_expected=sym(1);
    elseif mod(aud034_k,2),aud034_expected=sym(0);
    else,aud034_expected=-sym(1)/(aud034_k^2-1);end
    assert(isAlways(simplify(aud034_mean-aud034_expected)==0));
end
aud034_weight=0;
for aud034_i=0:127
    for aud034_j=0:aud034_i
        aud034_weight=aud034_weight+1+(aud034_i~=aud034_j);
    end
end
assert(aud034_weight==128^2);
aud034_eps=sym(9)/5000;
aud034_E=2^(-sym(4)/3)*sym(3)/4*((sym(256)/3)/sym(pi)^3)^(sym(1)/3) ...
    *((sym(16)/3)/sym(pi))^(sym(2)/3)*aud034_eps^(sym(4)/3);
aud034_allowed=((sym(29)/50000)^(sym(2)/3)-aud034_E)^(sym(3)/2);
assert(double(aud034_E)<.00013);
assert(double(aud034_allowed)>.0005638 && double(aud034_allowed)<.0005639);
fprintf('034 PASS exponents,orthant,chord,means,weighted_area1;Etime=%.17g allowed_IH=%.17g\n', ...
    double(aud034_E),double(aud034_allowed));
fprintf('034 theory only:full8256 I contribution certificate remains separate.\n');
