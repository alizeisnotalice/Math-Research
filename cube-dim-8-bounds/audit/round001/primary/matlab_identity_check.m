% Primary MATLAB R001 finite identity cross-check; no path or workspace reset.
cube8r001_rho=sym('cube8r001_rho','positive');
cube8r001_t=sym('cube8r001_t','positive');
cube8r001_q=1+1/(2*cube8r001_t);
cube8r001_antiderivative=cube8r001_q^7*log(cube8r001_rho);
for cube8r001_j=1:7
 cube8r001_antiderivative=cube8r001_antiderivative-nchoosek(7,cube8r001_j)*cube8r001_q^(7-cube8r001_j)*sym(-1/2)^cube8r001_j/(cube8r001_j*cube8r001_rho^cube8r001_j);
end
fprintf('n8 antiderivative identity: %d\n',isAlways(simplify(diff(cube8r001_antiderivative,cube8r001_rho)-(cube8r001_q-1/(2*cube8r001_rho))^7/cube8r001_rho)==0));
fprintf('short-time exponent = %s; integrated exponent = %s\n',char(sym(-6)/8+sym(7)/8),char(1+sym(-6)/8+sym(7)/8));
