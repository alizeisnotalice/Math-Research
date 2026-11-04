function round29_truncation()
% MATLAB independent Pareto-tail density integration; no workspace reset.
ns=[2,3]; gs=[.25,.5,.75]; cs=[.25,.5,.75,.95]; rs=[2,4,16,256]; maxErr=0; count=0;
for n=ns
    for g=gs
        for c=cs
            t=n*(1-g)*c; p=1/(1-t/n); ep=fullm(n,g,p); ell=ep^(1/p);
            for R=rs
                A=ell*R^(n-t); lep=fullm(n,g,p)-tailm(n,g,p,A); v=lep/ep+R^n*tailm(n,g,0,A);
                m=ell^(-t/(n-t))*lep+R^t*tailm(n,g,1,A); rel=m*v^(-t/n)/ell;
                assert(v>0 && v<=1 && rel>0 && rel<=1); count=count+1;
            end
            % Density-integral split audit for q=0,1,p at A=1,2,257.
            for q=[0,1,p]
                for A=[1,2,257]
                    dens=(n-1)/(n*g)*integral(@(x) x.^(q-1/g-1),1,A,'RelTol',1e-11,'AbsTol',1e-12);
                    maxErr=max(maxErr,abs(1/n+dens+tailm(n,g,q,A)-fullm(n,g,q)));
                end
            end
        end
    end
end
assert(maxErr<1e-9); fprintf('round29 MATLAB: %d finite-R rows; split error %.3g\n',count,maxErr);
end
function z=fullm(n,g,q), z=(n-g*q)/(n*(1-g*q)); end
function z=tailm(n,g,q,A)
e=q-1/g; z=(n-1)/(n*g)*A^e/(-e)*integral(@(y) exp(-y),0,Inf,'RelTol',1e-12,'AbsTol',1e-12);
end
