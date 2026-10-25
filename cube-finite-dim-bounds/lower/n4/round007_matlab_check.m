% Primary-agent serialized check; no clear/path/cwd changes.
n4r7sigma = .001; n4r7share = [.8 .2]; n4r7count = [16 24 18];
n4r7weights = [1 n4r7sigma*n4r7share];
n4r7mass = sum(n4r7weights.*n4r7count);
n4r7radius = n4r7mass^(1/4)/2;
n4r7S = sum(n4r7share.*n4r7count(2:3));
n4r7r0 = n4r7count(1)^(1/4)/2;
n4r7velocity = n4r7r0*n4r7S/(4*n4r7count(1));
assert(abs(n4r7mass/(2*n4r7radius)^4-1)<1e-12);
fprintf('n4 round007 joint mass=%.16g radius=%.16g cellmass=%.16g S=%.16g r1=%.16g\n',n4r7mass,n4r7radius,sum(n4r7weights),n4r7S,n4r7velocity);
