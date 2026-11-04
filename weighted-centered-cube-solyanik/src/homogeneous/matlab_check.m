% Independent adaptive integration of cube face formulas; numerical evidence.
checkRoot=fileparts(mfilename('fullpath'));
records=zeros(30,4); row=0;
for n=[2,3]
    for p=[1,2,4]
        for beta=[1/128,1/4,1/2,3/4,127/128]
            tau=n*beta;
            exponent=(tau-n)/p;
            if n==2
                I=integral(@(y)(1+y.^p).^exponent,0,1, ...
                    'AbsTol',1e-12,'RelTol',1e-12);
            else
                I=integral2(@(y,z)(1+y.^p+z.^p).^exponent, ...
                    0,1,0,1,'AbsTol',1e-12,'RelTol',1e-12);
            end
            Vp=2^n*gamma(1+1/p)^n/gamma(1+n/p);
            kappa=n*2^n/tau*I;
            S=(n/tau)*Vp^(1-beta)*(2^n)^beta/kappa;
            assert(isfinite(S)&&S>1);
            row=row+1;records(row,:)=[n,p,tau,S];
        end
    end
end
writematrix(records,fullfile(checkRoot,'matlab_selected.csv'));
fprintf('MATLAB %s: %d independent face-integral cases PASS.\n',version,row);
S=sqrt(pi)/(2*asinh(1));fprintf('n=2,p=2,tau=1: %.16g\n',S);
assert(abs(S-records(records(:,1)==2 & records(:,2)==2 & records(:,3)==1,4))<1e-11);
