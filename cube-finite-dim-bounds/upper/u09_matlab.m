% U09: one finite path + radial/normalization algebra only.
% Not a continuum kernel cost; preserve shared session (no clear/path).
upU09N=5;
upU09Q=diag(ones(upU09N-1,1),1)+diag(ones(upU09N-1,1),-1);
upU09Q=upU09Q-diag(sum(upU09Q,2));
upU09t=1/4; upU09P=eye(upU09N)+upU09t*upU09Q;
upU09G=[0;3/4;1/2;1/4;0]; upU09e=[0;1;0;0;0];
upU09c=3/8; upU09H=max(upU09G-upU09c,0);
upU09mu=upU09Q*upU09H+upU09e;
upU09slack=upU09P*upU09H-max(upU09P*upU09G-upU09c,0);
assert(norm(upU09Q-upU09Q',inf)==0);
assert(norm(sum(upU09Q,2),inf)==0);
assert(all(upU09P(:)>=0));
assert(all(upU09mu>=0) && abs(sum(upU09mu)-1)<1e-14);
assert(all(upU09slack>=0));
fprintf('U09 rows %.17g source_mass %.17g min_source %.17g\n',norm(sum(upU09P,2)-1,inf),sum(upU09mu),min(upU09mu));
fprintf('U09 raw_source');fprintf(' %.17g',upU09Q*upU09G);fprintf('\n');
fprintf('U09 positive_source');fprintf(' %.17g',upU09mu);fprintf('\n');
fprintf('U09 Jensen_slack');fprintf(' %.17g',upU09slack);fprintf('\n');
for upU09n=[3,4]
    upU09a=2;upU09cRad=3/4;
    upU09rho=(upU09a/upU09cRad)^(1/(upU09n-1));
    upU09mass=integral(@(upU09r)upU09a-upU09cRad*upU09r.^(upU09n-1),0,upU09rho);
    upU09expected=(upU09n-1)/upU09n*upU09cRad^(-1/(upU09n-1))*upU09a^(upU09n/(upU09n-1));
    % Synthetic g=J=1 checks only formula cancellation; no G-cost meaning.
    upU09cOpt=1/upU09n;
    upU09ratio=((upU09n-1)/upU09n*upU09cOpt^(-1/(upU09n-1)))/(2^upU09n*(1-upU09cOpt));
    upU09ratioExpected=upU09n^(1/(upU09n-1))/2^upU09n;
    assert(abs(upU09mass-upU09expected)<1e-12);
    assert(abs(upU09ratio-upU09ratioExpected)<1e-12);
    fprintf('U09 algebra n%d radial_defect %.17g normalized_formula_defect %.17g\n',upU09n,abs(upU09mass-upU09expected),abs(upU09ratio-upU09ratioExpected));
end
