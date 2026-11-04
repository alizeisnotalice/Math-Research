% Side-window check from the original four contact-scale values (double precision).
base=fileparts(mfilename('fullpath'));
ps=[-1.25 -1 -0.875 -0.75 -0.5 -0.25 -0.125 0 0.5 1];
ts=[0.02 0.01 0.003];
rows=zeros(2*numel(ps)*numel(ts),5);k=0;
for mode=1:2
    for p=ps
        for t=ts
            k=k+1;
            if mode==1
                beta=0.5+p*t;
                f=-4-16*p+G(1+4*p);
            else
                beta=1.5+p*t;
                f=12+16*p+G(-3-4*p);
            end
            diagonal=max(-expm1(log1p(-t)/beta),expm1(-log1p(-t)/(2-beta)))/2;
            waypoints=branchSwitches(mode,t,beta,diagonal/t);
            area=8*t^2*integral(@(s) arrayfun(@(u) radialWidth(u,t,beta),s), ...
                0,diagonal/t,'Waypoints',waypoints,'AbsTol',2e-13,'RelTol',2e-12);
            rows(k,:)=[mode p t area (area-4*t^2-f*t^3)/t^4];
        end
    end
end
assert(all(isfinite(rows),'all'));
assert(max(abs(rows(:,5)))<150);
result=struct('kind','independent-four-scale-boundary-integral', ...
    'version',version,'random_seed',[],'integral_AbsTol',2e-13, ...
    'integral_RelTol',2e-12,'fzero_TolX',2e-14, ...
    'waypoints','contact-value crossing roots from 65 scan nodes; not a root-count certificate', ...
    'columns',{{'mode_1_low_2_high','lambda','t','area','error_over_t4'}},'rows',rows);
fid=fopen(fullfile(base,'..','results','uniform_square_matlab.json'),'w');
assert(fid~=-1);fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('Side windows: %d independent four-scale boundary integrals.\n',k);
fprintf('Max abs error/t^4: %.12g\n',max(abs(rows(:,5))));

function value=G(q)
    if q<=0
        value=0;
    elseif q>=1
        value=2*q-2/3;
    else
        value=2*q^2-2*q^3/3;
    end
end

function width=radialWidth(s,t,beta)
    fun=@(r) max([beta*log1p(-2*t*r), ...
        (beta-1)*log1p(-2*t*s)+log1p(-t*(s+r)), ...
        (beta-2)*log1p(2*t*s)+log1p(t*(s-r)), ...
        (beta-2)*log1p(2*t*r)])-log1p(-t);
    if fun(s)<=0
        width=0;
    else
        width=fzero(fun,[s,0.49/t],optimset('TolX',2e-14,'Display','off'))-s;
    end
end

function roots=branchSwitches(mode,t,beta,top)
    if mode==1
        r=-expm1(log1p(-t)/beta)/(2*t);
        fun=@(s) (beta-1)*log1p(-2*t*s)+log1p(-t*(s+r))-log1p(-t);
    else
        r=expm1(-log1p(-t)/(2-beta))/(2*t);
        fun=@(s) (beta-2)*log1p(2*t*s)+log1p(t*(s-r))-log1p(-t);
    end
    nodes=linspace(0,top*(1-1e-8),65);
    values=arrayfun(fun,nodes);
    roots=zeros(1,numel(nodes));count=0;
    for j=1:numel(nodes)-1
        if values(j)*values(j+1)<0
            root=fzero(fun,nodes(j:j+1),optimset('TolX',2e-14,'Display','off'));
            count=count+1;roots(count)=root;
        end
    end
    roots=roots(1:count);
end
