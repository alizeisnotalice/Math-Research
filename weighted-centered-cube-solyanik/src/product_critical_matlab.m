% Numerical illustration only; theorem is proved in product_critical_extension.tex.
% w(x)=exp(x_1)*prod_{j>=2}(1+x_j^2)^(-1).
% T_1=exp(x_1)-1. F: T_1 in [0,1/m), other coordinates in [-1,1].
numChecksTopic2 = 0;
minGapTopic2 = Inf;
maxResidualTopic2 = 0;
for dimTopic2 = 1:6
 for arityTopic2 = [2 3 4 8]
  for integerMassTopic2 = [17 97]
   radiusTopic2 = asinh(integerMassTopic2/(2*exp(0.3)));
   assert(radiusTopic2 > 1.3);
   marginalTopic2 = atan(0.3+radiusTopic2)-atan(0.3-radiusTopic2);
   cubeMassTopic2 = integerMassTopic2*marginalTopic2^(dimTopic2-1);
   sharedMassTopic2 = (1/arityTopic2)*(pi/2)^(dimTopic2-1);
   baseTopic2 = 1-1/arityTopic2;
   gapTopic2 = baseTopic2^(arityTopic2-1)*sharedMassTopic2/cubeMassTopic2;
   assert(gapTopic2 > 0);
   residualTopic2 = abs(2*exp(0.3)*sinh(radiusTopic2)-integerMassTopic2);
   assert(residualTopic2 < 1e-10);
   minGapTopic2 = min(minGapTopic2,gapTopic2);
   maxResidualTopic2 = max(maxResidualTopic2,residualTopic2);
   numChecksTopic2 = numChecksTopic2+1;
  end
 end
end
fprintf('MATLAB %s\n',version);
fprintf('product weight samples=%d, n=1..6, m=2,3,4,8, N=17,97\n',numChecksTopic2);
fprintf('minimum positive gap=%.17g, maximum mass residual=%.17g\n',minGapTopic2,maxResidualTopic2);
fprintf('PASS: numerical samples only; exact integer-period argument is in proof.\n');
