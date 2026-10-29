% Independent Symbolic Math Toolbox check; x = epsilon^(1/n).
syms x
for n=[1 2 3 10]
    N=max(2*n+1,n+3);
    D=1/(1-(2*x/(1+2*x))^n);
    fprintf('n=%d d: %s\n',n,char(taylor(D,x,'Order',N+1)));
    for m=[1 2 3 10]
        % Expand the unit factor first, avoiding branch-dependent fractional powers.
        sm=m*(1-(1-x^n)^(sym(1)/m));
        unit=taylor(sm/x^n,x,'Order',N+1);
        y=x*taylor(unit^(sym(1)/n),x,'Order',N+1);
        c=taylor(subs(D,x,y),x,'Order',N+1);
        fprintf('n=%d m=%d: %s\n',n,m,char(c));
    end
end
