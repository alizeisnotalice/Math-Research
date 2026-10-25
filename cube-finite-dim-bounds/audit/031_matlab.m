% Parent serialized exact Gram and conditional target algebra only.
% No local functions, clear/path changes, Green evaluation or search.
for aud031_n=[3,4]
    aud031_body=sym([ones(2^(aud031_n-1),1),2*(dec2bin(0:2^(aud031_n-1)-1,aud031_n-1)-'0')-1]);
    aud031_face=sym(zeros(0,aud031_n));
    for aud031_i=1:aud031_n
        for aud031_j=aud031_i+1:aud031_n
            for aud031_sign=[-1,1]
                aud031_v=sym(zeros(1,aud031_n));
                aud031_v(aud031_i)=1;aud031_v(aud031_j)=aud031_sign;
                aud031_face=[aud031_face;aud031_v]; %#ok<AGROW>
            end
        end
    end
    aud031_mb=2^(aud031_n-1);aud031_mf=aud031_n*(aud031_n-1);
    assert(isequal(aud031_body.'*aud031_body,aud031_mb*sym(eye(aud031_n))));
    assert(isequal(aud031_face.'*aud031_face,2*(aud031_n-1)*sym(eye(aud031_n))));
    aud031_U=2^aud031_n;
    aud031_threshold=sym(aud031_U)^(aud031_n-1)/sym(aud031_n)^aud031_n;
    if aud031_n==3,assert(isAlways(aud031_threshold==sym(64)/27));
    else,assert(isAlways(aud031_threshold==16));end
    fprintf('031 Gram PASS n=%d body=%d face=%d\n',aud031_n,aud031_mb,aud031_mf);
end
aud031_g=sym(13)/2500;aud031_I=sym(57)/100000;
aud031_slack=64*aud031_g^3-27*aud031_I^2;
aud031_cost2=27*aud031_I^2/aud031_g^3;
assert(isAlways(aud031_slack==sym(56653)/250000000000));
assert(isAlways(aud031_cost2==sym(2193075)/35152));
assert(isAlways(aud031_cost2<64));
fprintf('031 conditional target PASS slack=%s cost2=%s cost=%s\n', ...
    char(aud031_slack),char(aud031_cost2),char(vpa(sqrt(aud031_cost2),50)));
fprintf('031 targets only: no g/I full-domain certificate and no new upper bound.\n');
