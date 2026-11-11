% r084 stable MATLAB floating sanity checks; batch execution exits on completion.
% No endpoint theorem or directed certificate depends on this script.
Ns=[1024 4096 16384 65536 1048576];
for N=Ns
    logJ=log(2)+(N-1)*log(1+log(2));
    logU=.75+(N/2)*log(pi)-gammaln(1+N/2)+(N/2)*log(6/N);
    L=-N*log(2)-2*log(N)-logU-N*log(8)-logJ;
    fprintf('N=%d log coarse lower=%0.9f\n',N,L);
end
fprintf('log2-2/3=%0.15g, sqrt(12*pi*exp(1))=%0.15g\n',log(2)-2/3,sqrt(12*pi*exp(1)));
