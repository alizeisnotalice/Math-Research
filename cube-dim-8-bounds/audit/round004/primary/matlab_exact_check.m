% Primary R004 MATLAB exact BB/peak cross-check, no workspace/path reset.
cube8r004_t=sym('cube8r004_t','positive');
cube8r004_PB=(25344*cube8r004_t^4+250*cube8r004_t^2+1)/(7077888*cube8r004_t^8+860160*cube8r004_t^6+33984*cube8r004_t^4+540*cube8r004_t^2+3);
cube8r004_peak=limit(cube8r004_t^4*cube8r004_PB,cube8r004_t,inf);
fprintf('P4 exact peak coefficient = %s\n',char(cube8r004_peak));
cube8r004_BB=int(cube8r004_PB^2,cube8r004_t,0,inf);
fprintf('BB integral / pi = %s\n',char(simplify(cube8r004_BB/sym(pi))));
fprintf('BB coefficient check = %d\n',isAlways(simplify(cube8r004_BB-sym(pi)*sym(143)/10584)==0));
