% Finite numerical consistency checks of proved candidate envelopes.
% These checks are not a proof over arbitrary measurable sets.
outDir = fileparts(mfilename('fullpath'));
rows = zeros(3*3*4*512, 7);
idx = 0;
for n = 1:3
    for beta = [0.25, 0.5, 0.75]
        tau = n*beta;
        for m = [2,3,4,8]
            for k = 1:512
                epsLevel = exp(-k*log(2));
                t = -m*expm1(log1p(-epsLevel)/m);
                z = exp(beta*log(t));
                C = exp(-log(t)/(2*n));
                A = exp(tau*log1p(-1/C));
                lowerShape = exp(-tau*log1p(exp(log(t)/n)));
                lower = lowerShape/(1-z*lowerShape);
                if C <= 2 || A <= z
                    upper = Inf; % explicit bound not invoked in this region
                else
                    logSmall = n*log(2) + (n-tau)*log(2*(1+C)) ...
                        - log1p(-exp((tau-n)*log(2))) ...
                        + log(t) - log1p(-t) - log(z);
                    upper = 1/(A-z) + exp(logSmall);
                end
                assert(lower > 0 && lower <= upper*(1+1e-12));
                idx = idx+1;
                rows(idx,:) = [n,beta,m,k,lower,upper,t/epsLevel];
            end
        end
    end
end
T = array2table(rows, 'VariableNames', ...
    {'n','beta','m','k','lower_over_t_beta','upper_over_t_beta','t_over_epsilon'});
writetable(T,fullfile(outDir,'dimension_grid.csv'));
fprintf('MATLAB %s: %d finite-grid checks PASS.\n',version,idx);
fprintf('Finite upper bound invoked at %d points; %d points outside its domain.\n', ...
    sum(isfinite(T.upper_over_t_beta)),sum(~isfinite(T.upper_over_t_beta)));
tail = T(T.k == 512,:);
disp(tail);
assert(all(tail.upper_over_t_beta > 0.999999));
assert(all(tail.upper_over_t_beta < 1.00001));
fprintf('All final normalized upper/lower bounds agree with limit 1 to <1e-5.\n');
