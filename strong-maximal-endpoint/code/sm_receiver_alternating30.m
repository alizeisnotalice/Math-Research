function result = sm_receiver_alternating30(n,N,K,outerIters,out)
%SM_RECEIVER_ALTERNATING30 Alternating finite-grid probe with independent ratios.
%
% This is a numerical initializer only, not an interval certificate.  The
% receiver centres are the r029 curved 2-block centres, but every point has
% its own (c_1,c_2).  Each outer sweep first ascends the *actual discretised*
% receiver average with H held fixed (using the exact kernel derivative), then
% solves the capped convex inner dual.  A proposed inner profile is accepted
% only when its strict discrete feasibility holds and its actual quadrature
% Orlicz cost is no larger than the retained old cost up to a 2e-10
% relative floating acceptance tolerance.
%
% Run, for example:
%   sm_receiver_alternating30(64,64,8,2)
%
% No conclusion about the continuum strong maximal operator follows from this
% script.  The dual profile is represented by positive kernel weights and is
% coordinate-nondecreasing formally, so exact strict point constraints would
% generate upper-orthant halo cones.  Here, however, those averages and the
% modular cost are only Gauss--Legendre estimates, not certified integrals.

if nargin < 1 || isempty(n), n = 64; end
if nargin < 2 || isempty(N), N = 64; end
if nargin < 3 || isempty(K), K = 8; end             % gives 1+2K receivers
if nargin < 4 || isempty(outerIters), outerIters = 2; end
if nargin < 5 || isempty(out)
    out = fullfile(tempdir,...
        sprintf('sm_receiver_alternating30_n%d_N%d_K%d.json',n,N,K));
end
assert(mod(n,2)==0 && n >= 8 && N >= 32 && K >= 1);

root = fileparts(fileparts(mfilename('fullpath')));
old = jsondecode(fileread(fullfile(root,'rounds','r026','checks.json')));
ref = jsondecode(fileread(fullfile(root,'rounds','r027','checks.json')));
base = old.rows([old.rows.n] == n & [old.rows.m] == 1);
assert(isscalar(base),'No r026 scalar baseline for this n.');
c0 = base.c;
rr = ref.rows;
isEqualTwo = reshape(arrayfun(@(r) numel(r.blocks)==2 && r.blocks(1)==r.blocks(2),rr),1,[]);
ix = find([rr.n] == n & isEqualTwo,1);
assert(~isempty(ix),'No r027 equal-block reference for this n.');
tauRef = rr(ix).tau;

b = n/2; A = 2; target = 1 + 1e-5;
mu = fzero(@(x) x - 1 - exp(-x),[1 2]);
spacing = 2*sqrt(n)/K;
v = (0:K)'*spacing;
d = .25*v.^2/(mu*n);
u = [0 0];
for j=2:K+1
    u = [u; d(j)+v(j), d(j)-v(j); d(j)-v(j), d(j)+v(j)]; %#ok<AGROW>
end
M = size(u,1);
cs = c0*ones(M,2);

% One common positive Gauss--Legendre grid.  It defines both the modular
% cost and all receiver averages, so retaining an old grid profile is exact
% within this numerical model.
p0 = 1/(1+exp(c0)); mz = 1+c0*p0; vz = 1+c0^2*p0*(1-p0);
lo = min(u,[],'all') - 1;
hi = max(u,[],'all') + b*mz + 16*sqrt(b*vz) + 1;
[t,qw] = localQuad(N,lo,hi);
[t1,t2] = ndgrid(t,t); [w1,w2] = ndgrid(qw,qw);
x1 = t1(:); x2 = t2(:); logQ = log(w1(:))+log(w2(:));
sumx = x1+x2; logCostWeight = logQ-sumx;
fprintf('r030 n=%d grid=%dx%d M=%d K=%d c0=%.12g\n',n,numel(t),numel(t),M,K,c0);

tic;
[R,P,mass] = kernelTables(cs);
if max(abs(mass-1)) >= 2e-6
    error('Quadrature range/mass inadequate: max error %.4g; masses %s.',...
        max(abs(mass-1)),mat2str(mass',5));
end
[H,w,st,inner] = solveOrRetain([],[],R,P,[]);
assert(inner.strict_scale_ok && min(P'*H) >= target-2e-10, ...
    'Initial grid profile is not strictly feasible.');
profileSourceCs = cs;
history = repmat(struct('iteration',0,'accepted_receiver_steps',0,...
    'max_receiver_average_gain',0,'old_min_average',0,'new_min_average',0,...
    'old_logJ',0,'new_logJ',0,'accepted_inner',false,'dual_exitflag',0,...
    'outer_grad_norm_inf',0,'ratio_min',0,'ratio_max',0),outerIters+1,1);
history(1) = makeHistory(0,0,0,NaN,st,st,true,inner,0,cs);

for sweep=1:outerIters
    Hold = H; stOld = st; Jold = st.J; minOld = min(P'*Hold);
    accepted = 0; maxGain = 0; gradInf = 0;
    % Each receiver has a distinct rectangle ratio.  Updating c_j only
    % changes its own average, so previous strict inequalities survive.
    for j=1:M
        [g,Bold] = receiverGradient(Hold,cs(j,:),u(j,:));
        gradInf = max(gradInf,max(abs(g)));
        ng = norm(g);
        if ng < 1e-12, continue; end
        dir = g/ng;
        acceptedHere = false;
        for step = [.30 .15 .075 .03 .01 .003]
            trial = min(3,max(.10,cs(j,:)+step*dir));
            Btrial = receiverAverage(Hold,trial,u(j,:));
            if Btrial >= Bold + 1e-10
                cs(j,:) = trial; accepted = accepted + 1;
                maxGain = max(maxGain,Btrial-Bold);
                acceptedHere = true;
                break;
            end
        end
        if ~acceptedHere
            % The exact derivative is only used to propose a direction; this
            % explicit recomputation is the acceptance test.
        end
    end
    [Rnew,Pnew,mass] = kernelTables(cs);
    assert(max(abs(mass-1)) < 2e-6);
    % Old H remains feasible after all accepted receiverwise ascents.
    oldUnderNew = Pnew'*Hold;
    if min(oldUnderNew) < minOld - 2e-9
        error('Outer acceptance invariant failed: old H lost feasibility.');
    end
    [Hcand,wcand,stcand,inner] = solveOrRetain(Hold,Jold,Rnew,Pnew,w);
    acceptedInner = inner.selected_new;
    if acceptedInner
        H=Hcand; w=wcand; st=stcand;
        profileSourceCs=cs;
    else
        % This is a genuine retained feasible discrete profile.  It protects
        % monotonicity against a local/failed fmincon solve.
        H=Hold; st=statsOf(H,Pnew); st.w = w;
    end
    P=Pnew;
    history(sweep+1) = makeHistory(sweep,accepted,maxGain,minOld,stOld,st,...
        acceptedInner,inner,gradInf,cs);
    fprintf(['sweep %d: receiver steps %d, max dB %.3e, minAvg %.9g, ',...
        'logJ %.10f, acceptedInner %d, c=[%.5g,%.5g]\n'],...
        sweep,accepted,maxGain,min(P'*H),st.logJ,acceptedInner,min(cs,[],'all'),max(cs,[],'all'));
end

union = discreteCurveUnion(u);
result = struct('matlab',version,'status',...
    ['finite positive quadrature / receiver-point alternating probe only; ',...
     'positive source weights formally define a monotone profile, but its ',...
     'continuum averages and modular cost are not interval-certified'],...
    'n',n,'block_size',b,'nodes_per_axis',numel(t),'quadrature_panels',ceil(N/16),...
    'receivers',M,'outer_sweeps',outerIters,'height_cap',A,'target_average',target,...
    'baseline_c',c0,'tau_ref',tauRef,'receiver_u',u,'receiver_c',cs,...
    'profile_source_c',profileSourceCs,...
    'quadrature_mass_max_error',max(abs(sum(P,1)'-1)),...
    'final_min_average',min(P'*H),'final_logJ',st.logJ,'final_J',st.J,...
    'curve_union_volume_diagnostic',union,'curve_union_log_score_diagnostic',log(union)-st.logJ,...
    'scalar_baseline_log_score',base.log_score,...
    'curve_minus_scalar_diagnostic',log(union)-st.logJ-base.log_score,...
    'profile_source_weights',w,'history',history,'seconds',toc);
fid=fopen(out,'w'); assert(fid>=0,'Cannot open output path.');
cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('FINAL score diagnostic %.10f; scalar %.10f; difference %+.8f; %.2fs\n',...
    result.curve_union_log_score_diagnostic,base.log_score,result.curve_minus_scalar_diagnostic,toc);

    function [Rloc,Ploc,masses] = kernelTables(cloc)
        Rloc=zeros(numel(sumx),M); Ploc=zeros(numel(sumx),M);
        for jj=1:M
            lr=logDensity(b,cloc(jj,1),x1-u(jj,1))+...
               logDensity(b,cloc(jj,2),x2-u(jj,2));
            Rloc(:,jj)=exp(lr+sum(u(jj,:))-tauRef);
            % logDensity is log(e^t f_{b,c}(t)); restore the two exponential
            % factors to obtain the probability-density quadrature weights.
            Ploc(:,jj)=exp(logQ+lr+sum(u(jj,:))-sumx);
        end
        masses=sum(Ploc,1)';
    end

    function [Hout,wout,stout,info] = solveOrRetain(Hold,Jhold,Rloc,Ploc,wstart)
        if nargin<5 || isempty(wstart), wstart=max(1e-4,ones(M,1)/M); end
        opts=optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true,...
            'Display','off','MaxIterations',350,'MaxFunctionEvaluations',6000,...
            'OptimalityTolerance',1e-8,'StepTolerance',1e-10);
        [wo,~,flag,solver]=fmincon(@dualObjective,wstart,[],[],[],[],...
            zeros(M,1),100*ones(M,1),[],opts);
        [Htry,wtry,stry,scaleOK]=makeStrict(wo,Rloc,Ploc);
        useNew=isempty(Hold) || (scaleOK && stry.J <= Jhold*(1+2e-10));
        info=struct('exitflag',flag,'iterations',solver.iterations,...
            'funcCount',solver.funcCount,'strict_scale_ok',scaleOK,...
            'selected_new',useNew);
        if useNew
            Hout=Htry; wout=wtry; stout=stry;
        else
            Hout=Hold; wout=wstart; stout=statsOf(Hold,Ploc);
        end

        function [f,g] = dualObjective(ww)
            a=Rloc*ww;
            hh=localProfile(log(max(a,realmin)),n,A);
            ee=Ploc'*hh;
            Jscaled=sum(exp(logCostWeight+tauRef+localLogPhi(hh,n)));
            dual=sum(ww)-ww'*ee+Jscaled;
            f=-dual; g=ee-1;
        end
    end

    function [hh,ww,ss,ok] = makeStrict(wo,Rloc,Ploc)
        loS=0; hiS=1; ok=true;
        for z=1:30
            hh0=localProfile(log(max(Rloc*(exp(hiS)*wo),realmin)),n,A);
            if min(Ploc'*hh0)>=target, break; end
            hiS=2*hiS;
            if hiS>128, ok=false; break; end
        end
        if ~ok, hh=zeros(size(x1)); ww=wo; ss=statsOf(hh,Ploc); return; end
        for z=1:45
            mid=(loS+hiS)/2;
            hh0=localProfile(log(max(Rloc*(exp(mid)*wo),realmin)),n,A);
            if min(Ploc'*hh0)<target, loS=mid; else, hiS=mid; end
        end
        ww=exp(hiS)*wo;
        hh=localProfile(log(max(Rloc*ww,realmin)),n,A);
        ss=statsOf(hh,Ploc);
    end

    function ss=statsOf(hh,Ploc)
        pp=Ploc'*hh;
        ll=localLogPhi(hh,n);
        JJ=sum(exp(logCostWeight+ll));
        ss=struct('averages',pp,'J',JJ,'logJ',log(JJ));
    end

    function [g,B] = receiverGradient(hh,cc,uu)
        B=receiverAverage(hh,cc,uu);
        pc=1./(1+exp(cc));
        b1=receiverAverageOrder(hh,cc(1),b-1,cc(2),b,uu(1)+cc(1),uu(2));
        b2=receiverAverageOrder(hh,cc(1),b,cc(2),b-1,uu(1),uu(2)+cc(2));
        g=b*pc.*[B-b1,B-b2];
    end

    function B=receiverAverage(hh,cc,uu)
        lr=logDensity(b,cc(1),x1-uu(1))+logDensity(b,cc(2),x2-uu(2));
        B=sum(exp(logQ+lr+sum(uu)-sumx).*hh);
    end

    function B=receiverAverageOrder(hh,cx,bx,cy,by,ux,uy)
        lr=logDensity(bx,cx,x1-ux)+logDensity(by,cy,x2-uy);
        B=sum(exp(logQ+lr+ux+uy-sumx).*hh);
    end
end

function h=localProfile(dd,n,A)
h=zeros(size(dd)); h(dd>=0)=1;
soft=dd>log(n);
y=max(expm1((dd(soft)-log(n))/(n-1)),0);
for j=1:8
    residual=(n-2)*log1p(y)+log(n+y)-dd(soft);
    deriv=(n-2)./(1+y)+1./(n+y);
    y=max(0,y-residual./deriv);
end
h(soft)=exp(min(y,log(A)));
end

function lp=localLogPhi(h,n)
lp=-Inf(size(h)); positive=h>0; y=log(h(positive));
lp(positive)=y+(n-1)*log1p(max(y,0));
end

function ld=logDensity(b,c,t)
% f_{b,c}(t)=sum_k binom(b,k)(t-ck)_+^(b-1)/(Gamma(b)(1+e^-c)^b).
if b==0
    error('The b=0 atom is not used here (n/2 >= 4).');
end
k=0:b; lw=gammaln(b+1)-gammaln(k+1)-gammaln(b-k+1);
xx=t-c*k; terms=lw+(b-1)*log(max(xx,realmin)); terms(xx<=0)=-Inf;
mx=max(terms,[],2); ld=mx+log(sum(exp(terms-mx),2))-gammaln(b)-b*log1p(exp(-c));
ld(t<=0)=-Inf;
end

function [t,w]=localQuad(N,lo,hi)
ii=(1:15)'; bb=ii./sqrt(4*ii.^2-1);
[V,D]=eig(diag(bb,1)+diag(bb,-1)); [x,ix]=sort(diag(D)); q=2*V(1,ix)'.^2;
panels=ceil(N/16); edges=linspace(lo,hi,panels+1);
half=diff(edges)/2; mid=(edges(1:end-1)+edges(2:end))/2;
t=x*half+mid; w=q*half; t=t(:); w=w(:);
end

function U=discreteCurveUnion(u)
% Diagnostic union of anchored cones with boundary points u.  It is not the
% physical witness volume when the receiver set is finite.
[a,ix]=sort(u(:,1)); b=cummin(u(ix,2));
U=sum(exp(-b(1:end-1)-a(1:end-1)).*(-expm1(-diff(a))))+exp(-a(end)-b(end));
end

function h=makeHistory(iter,steps,gain,minOld,stOld,stNew,acc,info,gInf,cs)
h=struct('iteration',iter,'accepted_receiver_steps',steps,'max_receiver_average_gain',gain,...
    'old_min_average',minOld,'new_min_average',min(stNew.averages),...
    'old_logJ',stOld.logJ,'new_logJ',stNew.logJ,'accepted_inner',acc,...
    'dual_exitflag',info.exitflag,'outer_grad_norm_inf',gInf,...
    'ratio_min',min(cs,[],'all'),'ratio_max',max(cs,[],'all'));
end
