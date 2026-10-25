frame8_V3=[1 -1 -1;1 -1 1;1 1 -1;1 1 1];
frame8_V4=[ones(8,1),repelem([-1;1],4),repmat(repelem([-1;1],2),2,1),repmat([-1;1],4,1)];
frame8_B3=[1 1 1;1 1 -1;1 -1 1]';
frame8_B4=[1 1 1 1;1 1 -1 -1;1 -1 1 -1;1 -1 -1 1]';
assert(isequal(frame8_V3'*frame8_V3,4*eye(3)));
assert(isequal(frame8_V4'*frame8_V4,8*eye(4)));
assert(abs(det(frame8_B3))==4 && abs(det(frame8_B4))==16);
fprintf('FRAME8 Gram3=4I Gram4=8I det3=%.17g det4=%.17g PASS\n',abs(det(frame8_B3)),abs(det(frame8_B4)));
