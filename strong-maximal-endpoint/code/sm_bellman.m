function result=sm_bellman(step,maxN,outputPath)
% Adaptive hit probability: exact Exp convolution of a linearly interpolated CDF.
% State interpolation is numerical; this is not a certified process bound.
if nargin<1, step=0.05; end
if nargin<2, maxN=256; end
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<3, outputPath=fullfile(root,'rounds','r011','bellman.json'); end
B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B); sigma=sqrt(mu);
q=0:step:ceil(maxN*mu+12*sqrt(maxN));
failure=tanh(q/2); % Exact one-step optimum: r=exp(-q).
decay=exp(-step); a=-expm1(-step); b=1-a/step;
rows=cell(0,1); snapshots=[2 4 8 16 32 64 128 256 512 1024 2048];
for n=2:maxN
    forcing=b*failure(2:end)+(a-b)*failure(1:end-1);
    conv=[0 filter(1,[1 -decay],forcing)];
    lower=zeros(size(q)); upper=q;
    for iter=1:38
        state=(lower+upper)/2; c=q-state; r=exp(-c);
        [convState,previousState]=evaluate(state,step,conv,failure);
        residual=conv+r.*convState-(1+r).*previousState;
        positive=residual>0;
        lower(positive)=state(positive); upper(~positive)=state(~positive);
    end
    state=(lower+upper)/2; r=exp(-(q-state));
    [convState,~]=evaluate(state,step,conv,failure);
    next=(conv+r.*convState)./(1+r); next(1)=0;
    assert(min(diff(next))>-1e-11,'Bellman monotonicity failed');
    assert(min(next)>-1e-12 && max(next)<1+1e-12);
    failure=next;
    if ismember(n,snapshots)
        beta=sigma*sqrt(log(n)/n);
        index=find(failure>=beta,1); assert(~isempty(index));
        adaptiveQ=q(index-1)+step*(beta-failure(index-1))/(failure(index)-failure(index-1));
        lo=0; hi=n*mu+12*sqrt(n);
        for k=1:70
            mid=(lo+hi)/2; [~,~,value]=sm_tail(n,1/B,mid);
            if value<beta, lo=mid; else, hi=mid; end
        end
        fixedQ=(lo+hi)/2;
        probe=n*mu-sigma*sqrt(n*log(n));
        [~,~,fixedFailure]=sm_tail(n,1/B,probe);
        adaptiveFailure=interp1(q,failure,probe,'linear');
        row=struct('n',n,'step',step,'soft_failure_target',beta, ...
            'adaptive_quantile',adaptiveQ,'fixed_rstar_quantile',fixedQ, ...
            'quantile_gain',adaptiveQ-fixedQ,'probe_q',probe, ...
            'adaptive_failure',adaptiveFailure,'fixed_failure',fixedFailure);
        rows{end+1}=row; %#ok<AGROW>
        fprintf('Bellman n=%3d h=%.4g quantileGain=%.9g beta=%.6g failRatio=%.6g\n', ...
            n,step,row.quantile_gain,beta,adaptiveFailure/fixedFailure);
    end
end
result=struct('matlab_version',version,'method', ...
    'piecewise-linear state interpolation; exact exponential cell convolution; unique-control bisection', ...
    'precision','double; no directed error certificate','step',step,'q_max',q(end),'rows',{rows});
fid=fopen(outputPath,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end

function [convState,previousState]=evaluate(s,step,conv,previous)
j=min(floor(s/step)+1,numel(previous)-1); h=s-(j-1)*step;
slope=(previous(j+1)-previous(j))/step;
a=-expm1(-h);
convState=exp(-h).*conv(j)+a.*previous(j)+max(h-a,0).*slope;
previousState=previous(j)+h.*slope;
end
