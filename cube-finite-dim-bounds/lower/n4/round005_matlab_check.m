% Shared-session-safe local velocity-interval checks. No clear/path/cwd changes.
n4r5L = [-0.8, -0.2, 0.3, -0.5];
n4r5R = [0.5, 0.1, -0.3, 0.8];
n4r5B = [0.1, 0.4, 0, 0.2];
n4r5T = 4;
for n4r5i = 1:numel(n4r5L)
    n4r5gap = max(0,n4r5R(n4r5i)-n4r5L(n4r5i)) ...
        - max(0,min(n4r5R(n4r5i),n4r5B(n4r5i))-max(n4r5L(n4r5i),-n4r5B(n4r5i)));
    n4r5intervals = sortrows([-n4r5T,n4r5L(n4r5i); n4r5R(n4r5i),n4r5T; -n4r5B(n4r5i),n4r5B(n4r5i)]);
    n4r5lo = n4r5intervals(1,1); n4r5hi = n4r5intervals(1,2); n4r5covered = 0;
    for n4r5k = 2:3
        if n4r5intervals(n4r5k,1) > n4r5hi
            n4r5covered = n4r5covered + n4r5hi-n4r5lo;
            n4r5lo = n4r5intervals(n4r5k,1); n4r5hi = n4r5intervals(n4r5k,2);
        else
            n4r5hi = max(n4r5hi,n4r5intervals(n4r5k,2));
        end
    end
    n4r5covered = n4r5covered+n4r5hi-n4r5lo;
    assert(abs((2*n4r5T-n4r5covered)-n4r5gap)<1e-12);
    fprintf('n4 round005 seam case %d exact-form gap=%.16g interval-union gap=%.16g\n',n4r5i,n4r5gap,2*n4r5T-n4r5covered);
end
