% U03 independent selected killed Levy matrix, not fractional power of a Laplacian.
upU03alpha=1.8;
upU03N=129;
upU03h=2/upU03N;
upU03x=((0:upU03N-1)-(upU03N-1)/2)*upU03h;
upU03k=1:upU03N-1;
upU03rates=upU03h^(-upU03alpha)/upU03alpha*((upU03k-.5).^(-upU03alpha)-(upU03k+.5).^(-upU03alpha));
upU03Q=toeplitz([-2^(upU03alpha+1)*upU03h^(-upU03alpha)/upU03alpha,upU03rates]);
upU03killing=((1-upU03x).^(-upU03alpha)+(1+upU03x).^(-upU03alpha))/upU03alpha;
[upU03V,upU03L]=eig(-upU03Q,'vector');
[upU03lambda,upU03ord]=sort(upU03L);
upU03V=upU03V(:,upU03ord);
upU03source=(upU03N+1)/2;
upU03scaled=upU03lambda/upU03lambda(1);
upU03scoeff=upU03V(upU03source,:).*sum(upU03V,1);
fprintf('U03 alpha=%.1f N=%d row killing defect=%.5g lambda0=%.17g\n',upU03alpha,upU03N,max(abs(sum(upU03Q,2)+upU03killing')),upU03lambda(1));
for upU03n=[3,4]
    upU03offset=35-upU03n;
    upU03a=upU03offset*upU03h;
    upU03pcoeff=upU03V(upU03source,:).*upU03V(upU03source+upU03offset,:)/upU03h;
    upU03mass=integral(@(upU03t) arrayfun(@(upU03z) (upU03scoeff*exp(-upU03scaled*upU03z))^upU03n,upU03t),0,Inf,'AbsTol',1e-11,'RelTol',1e-10)/upU03lambda(1);
    upU03green=integral(@(upU03t) arrayfun(@(upU03z) (upU03pcoeff*exp(-upU03scaled*upU03z))^upU03n,upU03t),0,Inf,'AbsTol',1e-11,'RelTol',1e-10)/upU03lambda(1);
    upU03cost=upU03mass/(2^upU03n*upU03a^upU03n*upU03green);
    fprintf('U03 n=%d a=%.17g mass=%.17g green=%.17g cost=%.17g\n',upU03n,upU03a,upU03mass,upU03green,upU03cost);
end
