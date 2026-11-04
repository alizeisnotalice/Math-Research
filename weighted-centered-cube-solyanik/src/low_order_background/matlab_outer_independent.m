% Independent floating-point outer-halo check.  No clearvars/path reset: this
% is deliberately safe to run in the connected user's MATLAB session.
% H_alpha is determined solely by Lebesgue length, not by the weight.
r21_root = fileparts(mfilename('fullpath'));
r21_count = 0; r21_fail = 0; r21_maxSlack = -Inf;
r21_maxFinalDifference = 0; r21_examples = struct([]);
for r21_case = 1:10
    r21_z = [-0.60 - 0.04*mod(r21_case,3), 0.60 + 0.03*mod(r21_case,4)];
    r21_tau = [0.15 + 0.10*mod(r21_case,5), 0.25 + 0.10*mod(r21_case+2,5)];
    r21_c = [0.45 + 0.05*mod(r21_case,4), 0.60 + 0.06*mod(r21_case,3)];
    % b=1+sum h 1_(a,b), hence 1<=b<=1+max(h) everywhere.
    r21_steps = [-3 -1 1; -1 0 3+2*mod(r21_case,4); 0 1 1+mod(r21_case,3); 1 3 1];
    for r21_j = 1:2
        for r21_k = [5 9 12]
            r21_eps = 2^(-r21_k);
            r21_E = [r21_z(r21_j)-1.5*r21_eps, r21_z(r21_j)-0.05*r21_eps; ...
                      r21_z(r21_j)+0.05*r21_eps, r21_z(r21_j)+1.5*r21_eps];
            r21_H = r21_outer_halo(r21_E,0.90);
            r21_me = r21_masses(r21_E,r21_steps,r21_c,r21_z,r21_tau);
            r21_mh = r21_masses(r21_H,r21_steps,r21_c,r21_z,r21_tau);
            r21_ratios = r21_mh./r21_me;
            r21_mix = sum(r21_mh)/sum(r21_me);
            r21_slack = r21_mix-max(r21_ratios);
            r21_pure = r21_ratios(r21_j+1);
            r21_diff = r21_mix-r21_pure;
            r21_count = r21_count+1;
            r21_fail = r21_fail + (r21_slack > 2e-11);
            r21_maxSlack = max(r21_maxSlack,r21_slack);
            if r21_k==12, r21_maxFinalDifference=max(r21_maxFinalDifference,abs(r21_diff)); end
            r21_examples(r21_count).case = r21_case; %#ok<SAGROW>
            r21_examples(r21_count).site = r21_j;
            r21_examples(r21_count).k = r21_k;
            r21_examples(r21_count).mixture = r21_mix;
            r21_examples(r21_count).componentMaximum = max(r21_ratios);
            r21_examples(r21_count).convexSlack = r21_slack;
            r21_examples(r21_count).pureAtSite = r21_pure;
            r21_examples(r21_count).localDifference = r21_diff;
        end
    end
end
r21_summary = struct('parameter_version','matlab-independent-outer-v1-20261105', ...
    'examples',r21_count,'failures',r21_fail,'max_convex_slack',r21_maxSlack, ...
    'max_abs_local_difference_k12',r21_maxFinalDifference,'matlab_version',version, ...
    'records',{r21_examples});
r21_fid=fopen(fullfile(r21_root,'matlab_outer_summary.json'),'w');
assert(r21_fid>=0,'Cannot open MATLAB summary for writing.');
fprintf(r21_fid,'%s\n',jsonencode(r21_summary)); fclose(r21_fid);
assert(r21_fail==0,'A convex-combination check failed.');
fprintf('MATLAB outer independent: %d examples, failures=%d, max slack %.3g\n',r21_count,r21_fail,r21_maxSlack);

function H=r21_outer_halo(E,alpha)
% Independent floating-point counterpart of the endpoint/cell construction.
ends=reshape(E.',1,[]); ell=sum(E(:,2)-E(:,1));
lo=ends(1)-ell/(2*alpha); hi=ends(end)+ell/(2*alpha); P=E;
for e=ends
    cuts=unique([lo hi (e+ends)/2]); cuts=cuts(cuts>lo & cuts<hi); cuts=[lo cuts hi];
    for q=1:numel(cuts)-1
        a=cuts(q); b=cuts(q+1); ga=r21_G(a,e,E,alpha); gb=r21_G(b,e,E,alpha);
        if ga>0 && gb>0, P=[P;a b]; %#ok<AGROW>
        elseif ga>0, P=[P;a, a+(b-a)*ga/(ga-gb)]; %#ok<AGROW>
        elseif gb>0, P=[P;a+(b-a)*(-ga)/(gb-ga), b]; %#ok<AGROW>
        end
    end
end
P=sortrows(P,1); H=[];
for q=1:size(P,1)
    if isempty(H) || P(q,1)>H(end,2)+1e-13, H=[H;P(q,:)]; %#ok<AGROW>
    else, H(end,2)=max(H(end,2),P(q,2)); end
end
end

function g=r21_G(x,e,E,alpha)
l=min(e,2*x-e); r=max(e,2*x-e);
g=sum(max(0,min(r,E(:,2))-max(l,E(:,1))))-2*alpha*abs(x-e);
end

function masses=r21_masses(I,steps,c,z,tau)
base=sum(I(:,2)-I(:,1)); bg=base;
for q=1:size(steps,1)
    bg=bg+steps(q,3)*sum(max(0,min(I(:,2),steps(q,2))-max(I(:,1),steps(q,1))));
end
masses=[bg zeros(1,numel(c))];
for j=1:numel(c)
    masses(j+1)=c(j)/tau(j)*sum(r21_sp(I(:,2)-z(j),tau(j))-r21_sp(I(:,1)-z(j),tau(j)));
end
end

function y=r21_sp(x,p)
y=sign(x).*abs(x).^p;
end
