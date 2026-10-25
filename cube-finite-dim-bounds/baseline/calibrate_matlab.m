% Numerical cross-check only; exact proof is in calibration_and_transfer.md.
format long g
for n = [3,4]
    u = 2^(1/n); a = (1+u)/2;
    An = 1.5*a^(n-1)-0.5*((3*a-2*(1.5)^(1/n))^(n-1)-((3-u)/2)^(n-1));
    fprintf('n=%d Aldaz=%.16g\n',n,An);
end
fprintf('radial3=%.16g radial4=%.16g\n',27*pi*sqrt(3)/16,2*pi^2);
