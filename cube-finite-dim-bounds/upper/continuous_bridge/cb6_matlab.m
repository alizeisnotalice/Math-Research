% cb6 finite normalization/outer-margin cross-check, no clear/path changes.
cb6_nlist=[3,4];
cb6_ranges=[-4/15,2/5;-1/4,1/2];
for cb6_idx=1:2
    cb6_n=cb6_nlist(cb6_idx);
    cb6_q=2*(cb6_n-1)/(cb6_n*(cb6_n+2));
    cb6_p=(cb6_n-1)/(cb6_n+2);
    cb6_axis=[1,zeros(1,cb6_n-1)];
    cb6_diag=ones(1,cb6_n)/sqrt(cb6_n);
    cb6_Haxis=sum(cb6_axis.^4)-3*sum(cb6_axis.^2)^2/(cb6_n+2);
    cb6_Hdiag=sum(cb6_diag.^4)-3*sum(cb6_diag.^2)^2/(cb6_n+2);
    assert(abs(cb6_Haxis-cb6_p)<1e-14);
    assert(abs(cb6_Hdiag+cb6_q)<1e-14);
    cb6_outermax=2^(2-cb6_n)-1+cb6_q*2^4/100;
    cb6_sigma=2*pi^(cb6_n/2)/gamma(cb6_n/2);
    cb6_cn=1/((cb6_n-2)*cb6_sigma);
    cb6_radial_mass=cb6_cn*cb6_sigma*(1/2-1/cb6_n);
    cb6_source_mass=cb6_cn*(cb6_n-2)*cb6_sigma;
    assert(cb6_outermax<0);
    assert(abs(cb6_radial_mass-1/(2*cb6_n))<1e-14);
    assert(abs(cb6_source_mass-1)<1e-14);
    fprintf('cb6 n=%d PASS H4=[%.16g,%.16g], outer max %.16g, radial mass %.16g, source mass %.16g\n',...
      cb6_n,-cb6_q,cb6_p,cb6_outermax,cb6_radial_mass,cb6_source_mass);
end
% These are legality small examples; neither radial example dominates cube.
