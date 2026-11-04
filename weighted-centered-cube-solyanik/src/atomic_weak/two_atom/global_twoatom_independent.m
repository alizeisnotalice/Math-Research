function global_twoatom_independent()
% Independent double-precision diagnostic.  It uses u=(r-1)/(r+1), not F(r)
% bisection.  It is numerical evidence, not a substitute for the sign proof.
taus = [1e-12,1e-9,1e-6,1e-4,.001,.01,.05,.1,.25,.5,.75,.9,.99,.999999,1-1e-12];
out = struct('tau',cell(numel(taus),1),'u',[],'r',[],'tau_J',[],'tau_R',[], ...
    'tau_gap',[],'K_residual',[],'in_interval',[],'strict_gap',[],'status',[]);
for i = 1:numel(taus)
    t = taus(i); upper = t/(2-t);
    lo = upper*1e-13; hi = upper*(1-1e-13);
    klo = localK(lo,t); khi = localK(hi,t);
    if ~(isfinite(klo) && isfinite(khi) && klo > 0 && khi < 0)
        out(i) = localrow(t,NaN,NaN,NaN,NaN,NaN,NaN,false,false,'bracket_failure');
        continue
    end
    for k=1:180
        mid=(lo+hi)/2;
        if localK(mid,t)>0, lo=mid; else, hi=mid; end
    end
    u=(lo+hi)/2; r=(1+u)/(1-u); d=r-1;
    J=(d^t+(2+d)^t)/(2*t*(1+d));
    R=((1+2*d)^t-(1-2*d)^t+d^t+(2-d)^t)/(2*t*(1+d));
    % Pairwise expm1 reduces, but cannot completely remove, cancellation when t is tiny.
    g1=exp(t*log(1-2*d))*expm1(t*(log(1+2*d)-log(1-2*d)));
    g2=exp(t*log(2+d))*expm1(t*(log(2-d)-log(2+d)));
    gap=(g1+g2)/(2*(1+d)); % equals t*(R-J) algebraically
    resolvable = abs(gap) > 128*eps(max(abs(g1),abs(g2)));
    nearneutral = abs(1-t) < 1e-10;
    if resolvable && ~nearneutral
        status='resolved_double'; strict=gap>0;
    elseif resolvable
        status='gap_resolved_root_flat'; strict=gap>0;
    else
        status='indeterminate_double_cancellation'; strict=NaN;
    end
    out(i)=localrow(t,u,r,t*J,t*R,gap,localK(u,t),r>1 && r<4/3,strict,status);
end
payload=struct('version','global-twoatom-matlab-v1-20261105','method', ...
 'log-u-bisection; double-precision diagnostic only','rows',{out}, ...
 'matlab_version',version,'timestamp',char(datetime('now','Format',"yyyy-MM-dd'T'HH:mm:ss")));
base=fileparts(mfilename('fullpath'));
fid=fopen(fullfile(base,'global_twoatom_matlab_summary.json'),'w');
fwrite(fid,jsonencode(payload),'char'); fclose(fid);
resolved = strcmp({out.status},'resolved_double') | strcmp({out.status},'gap_resolved_root_flat');
if ~all([out.in_interval]) || ~all([out(resolved).strict_gap]), error('global_twoatom:failed','A MATLAB resolved diagnostic check failed.'); end
fprintf('global_twoatom MATLAB rows=%d, resolved=%d, indeterminate=%d; all resolved interval/gap checks passed\n',numel(out),sum(resolved),sum(~resolved));
end
function z=localK(u,t)
z=-(1-t)*log(u)+log(t-(2-t)*u)-log((2-t)-t*u);
end
function q=localrow(t,u,r,tj,tr,gap,res,inside,strict,status)
q=struct('tau',t,'u',u,'r',r,'tau_J',tj,'tau_R',tr,'tau_gap',gap, ...
 'K_residual',res,'in_interval',inside,'strict_gap',strict,'status',status);
end
