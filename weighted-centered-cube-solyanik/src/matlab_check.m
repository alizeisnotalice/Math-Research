% Independent MATLAB R2026a floating-point check for Topic 2 scalar data.
ms = [2 3 4 8];
ns = [1 2 3];
epsilons = [1e-9 1e-6 1e-3];
alphas = [0.5 0.75 0.9 0.99 0.999];
fprintf('MATLAB %s\n', version);
fprintf('critical-side checks\n');
for m = ms
    star = ((m-1)/m)^m;
    for e = epsilons
        lambdas = [star*(1-e), star+(1-star)*e];
        labels = {'below','above'};
        for k = 1:2
            s = m*(1-lambdas(k)^(1/m));
            phi = 1-s;
            expected = (-1)^(k); % below: negative, above: positive
            pass = sign(phi) == expected;
            fprintf('m=%d %s eps=%g phi=%+.16e pass=%d\n',...
                m, labels{k}, e, phi, pass);
        end
    end
end
fprintf('star ratios\n');
for n = ns
    for alpha = alphas
        target = 1 + 2^n*(1-alpha)/alpha;
        for eta = [1e-2 1e-4 1e-6]
            q = (1-eta)*(1-alpha)/alpha;
            ratio = 1 + 2^n*q;
            strictAverage = 1/(1+q);
            pass = (ratio < target) && (strictAverage > alpha);
            fprintf('n=%d alpha=%.12g eta=%g error=%.3e pass=%d\n',...
                n, alpha, eta, target-ratio, pass);
        end
    end
end
