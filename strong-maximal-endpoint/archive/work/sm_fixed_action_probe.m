function sm_fixed_action_probe
addpath('/Users/zhengzhihao/Desktop/T/Math-Research-strong-maximal/strong-maximal-endpoint/code');
B=fzero(@(b)log(b)-1-1/b,[3 4]); mu=log(B); sigma=sqrt(mu);
for n=[16 32 64 128 256 512 1024 2048 4096]
 beta=sigma*sqrt(log(n)/n);
 q0=quant(n,mu,beta,mu);
 [c,negq,flag]=fminbnd(@(c)-quant(n,c,beta,mu),.3,2.3,optimset('TolX',1e-8));
 fprintf('n=%4d c=%.10g qgain=%.10g gain/logn=%.8g cshift/rho=%.8g flag=%d\n',n,c,-negq-q0,(-negq-q0)/log(n),(c-mu)/sqrt(log(n)/n),flag);
end
end
function q=quant(n,c,beta,mu)
 lo=0; hi=n*mu+12*sqrt(n);
 for j=1:55
  mid=(lo+hi)/2; [~,~,failure]=sm_tail(n,exp(-c),mid);
  if failure<beta,lo=mid;else,hi=mid;end
 end
 q=(lo+hi)/2;
end
