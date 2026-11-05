function sm_checks
% Analytic distribution controls, not optimization self-comparisons.
for n = [1 2 8 64 256]
    for q = [0 0.1 n 2*n]
        [a,~,b] = sm_tail(n,1,q);
        assert(abs(a-gammainc(q,n,'upper'))<2e-13);
        assert(abs(a+b-1)<2e-13);
        assert(abs(sm_tail(n,0,q)-gammainc(q,n,'upper'))<2e-13);
    end
end
for r = [0.1 0.3 0.8]
    q = -log(r)+2;
    assert(abs(sm_tail(1,r,q)-2*exp(-q)/(1+r))<2e-14);
end
fprintf('PASS: r=0,1 Gamma controls; n=1 exact tail; complementary probabilities.\n');
end
