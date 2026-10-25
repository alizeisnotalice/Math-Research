% Conditional mass algebra cross-check, not actual Green kernel certification.
for mass9_n=[3,4]
 mass9_g=1.7; mass9_c=mass9_g/mass9_n;
 mass9_z=(mass9_g/mass9_c)^(1/(mass9_n-1));
 mass9_num=integral(@(mass9_s) mass9_g-mass9_c*mass9_s.^(mass9_n-1),0,mass9_z,'AbsTol',1e-12,'RelTol',1e-12);
 mass9_exact=(mass9_n-1)/mass9_n*mass9_g^(mass9_n/(mass9_n-1))*mass9_c^(-1/(mass9_n-1));
 assert(abs(mass9_num-mass9_exact)<1e-11);
 mass9_threshold=2^mass9_n/mass9_n^(mass9_n/(mass9_n-1));
 fprintf('MASS9 n=%d radial_error=%.17g required_average=%.17g PASS\n',mass9_n,abs(mass9_num-mass9_exact),mass9_threshold);
end
mass9_a=0.5;
mass9_axis=integral(@(mass9_t) mass9_t.*(mass9_a^2+20*mass9_t.^2)./(4*pi^3*(mass9_a^2+4*mass9_t.^2).^3),0,Inf,'AbsTol',1e-13,'RelTol',1e-12);
assert(abs(mass9_axis-3/(8*pi^3))<1e-12);
fprintf('MASS9 axis_G=%.17g error=%.17g PASS\n',mass9_axis,abs(mass9_axis-3/(8*pi^3)));
