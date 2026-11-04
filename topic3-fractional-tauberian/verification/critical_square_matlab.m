% Independent geometry: optimize the four exact contact scales, then integrate.
% Double precision is used; asymptotic residuals are numerical evidence only.
base=fileparts(mfilename('fullpath'));
ts=[0.02 0.01 0.003];
parameters={[-2 -1 -0.75 -0.5 -0.25 0 0.5 1 2], ...
    [-1 -0.6 -0.5 -0.45 -0.4 -0.375 -0.3 -0.25 -0.1 0 1]};
rows=zeros(sum(cellfun(@numel,parameters))*numel(ts),6);
k=0;
for mode=1:2
    for p=parameters{mode}
        for t=ts
            k=k+1;
            if mode==1
                beta=1+p*t;
                approximation=2*t^2+max(-2*p,2+2*p)*t^3;
                power=4;
            else
                beta=1-t/2+p*t^2;
                approximation=2*t^2+t^3+innerQ(p)*t^4;
                power=5;
            end
            diagonal=max(-expm1(log1p(-t)/beta),expm1(-log1p(-t)/(2-beta)))/2;
            top=diagonal/t;
            area=8*t^2*integral(@(s) arrayfun(@(u) radialWidth(u,t,beta),s), ...
                0,top,'AbsTol',2e-13,'RelTol',2e-12);
            rows(k,:)=[mode p t area (area-approximation)/t^power top];
        end
    end
end
assert(all(isfinite(rows),'all'));
assert(max(abs(rows(rows(:,1)==1,5)))<20);
assert(max(abs(rows(rows(:,1)==2,5)))<5);
result=struct('kind','independent-double-precision-numerical-evidence', ...
    'version',version,'random_seed',[],'integral_AbsTol',2e-13, ...
    'integral_RelTol',2e-12,'fzero_TolX',2e-14, ...
    'columns',{{'mode_1_outer_2_inner','parameter','t','area','scaled_window_error','scaled_diagonal'}}, ...
    'rows',rows);
fid=fopen(fullfile(base,'..','results','critical_square_matlab.json'),'w');
assert(fid~=-1);fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('Critical square: %d independent four-scale boundary integrals.\n',k);
fprintf('Outer scaled residual max %.9g, inner %.9g\n', ...
    max(abs(rows(rows(:,1)==1,5))),max(abs(rows(rows(:,1)==2,5))));

function width=radialWidth(s,t,beta)
    fun=@(r) max([beta*log1p(-2*t*r), ...
        (beta-1)*log1p(-2*t*s)+log1p(-t*(s+r)), ...
        (beta-2)*log1p(2*t*s)+log1p(t*(s-r)), ...
        (beta-2)*log1p(2*t*r)])-log1p(-t);
    if fun(s)<=0
        width=0;
    else
        root=fzero(fun,[s,0.49/t],optimset('TolX',2e-14,'Display','off'));
        width=root-s;
    end
end

function q=innerQ(mu)
    a=1+2*mu;
    if a<=0
        q=-1/6-2*mu;
    elseif a>=1/2
        q=7/6+2*mu;
    else
        q=5/6-a+8*a^3/3;
    end
end
