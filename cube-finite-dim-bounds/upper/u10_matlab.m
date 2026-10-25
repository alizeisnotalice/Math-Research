% U10 fixed support/mass/tail arithmetic only. No candidate, LP, G/mu integral.
% Preserve shared session; upU10 prefix; no clear/path.
upU10centers=[3/4*eye(3);-3/4*eye(3)];upU10d2=[];
for upU10i=1:6
    for upU10j=upU10i+1:6
        upU10d2(end+1)=sum((upU10centers(upU10i,:)-upU10centers(upU10j,:)).^2);
    end
end
assert(min(upU10d2)>1);
upU10M0=integral(@(upU10r)4*pi*(1-upU10r.^2).^2.*upU10r.^2,0,1);
upU10M1=6*(1/2)^3*upU10M0;
assert(abs(upU10M0-32*pi/105)<1e-12 && abs(upU10M1-8*pi/35)<1e-12);
upU10V=[1 1 1;1 1 -1;1 -1 1;1 -1 -1];
upU10B=upU10V(1:3,:)';
assert(abs(abs(det(upU10B))-4)<1e-14);
upU10C=1/(4*pi^3);upU10delta=1/(8*sqrt(3));
upU10A=6*upU10C/upU10delta^4;
upU10cross=(upU10C/upU10A)^(1/4);
upU10Gbound=upU10A*upU10cross^2/2+upU10C/(2*upU10cross^2);
assert(abs(upU10Gbound-48*sqrt(6)/pi^3)<1e-12);
upU10v=upU10V(1,:);upU10w=upU10V(2,:);upU10shifted=[];
for upU10i=1:4
    upU10a2=sum(upU10w.^2)-(upU10w*upU10V(upU10i,:)')^2/3;
    upU10a1=2*(upU10v*upU10w'-(upU10v*upU10V(upU10i,:)')*(upU10w*upU10V(upU10i,:)')/3);
    upU10a0=sum(upU10v.^2)-(upU10v*upU10V(upU10i,:)')^2/3;
    upU10shifted(end+1,:)=[upU10a2,2*upU10a2+upU10a1,upU10a2+upU10a1+upU10a0-8/3];
end
assert(min(upU10shifted(:))>-1e-12);
% Synthetic weights solely test the closed tail formula, not feasibility.
upU10c=(4*log(2)-1)/(24*pi^3);upU10W=2/1000;
upU10T=max(60,160*upU10W/(sqrt(3)*upU10c));
upU10R=sqrt(3)*upU10T+5/4;
upU10tailSlack=upU10c/(4*pi*upU10T)-40*upU10W/(pi*sqrt(3)*upU10T^2);
assert(upU10tailSlack>=-1e-14);
fprintf('U10 supports min_distance_sq %.17g outer_radius 1.25 Omega_buffer .25\n',min(upU10d2));
fprintf('U10 M0 %.17g M1 %.17g ratio %.17g mass_defect %.17g\n',upU10M0,upU10M1,upU10M1/upU10M0,abs(upU10M0-32*pi/105));
fprintf('U10 basis_absdet %.17g G_bound_at_d1 %.17g defect %.17g\n',abs(det(upU10B)),upU10Gbound,abs(upU10Gbound-48*sqrt(6)/pi^3));
fprintf('U10 geometry shifted_coeff_min %.17g\n',min(upU10shifted(:)));
fprintf('U10 synthetic_only c %.17g T %.17g R %.17g tail_slack %.17g; not a candidate\n',upU10c,upU10T,upU10R,upU10tailSlack);
