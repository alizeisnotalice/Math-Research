% Fixed analytic interface cross-check, no corrected-kernel feasibility claim.
for bump9_n=[3,4]
 bump9_bits=dec2bin(0:2^(bump9_n-1)-1)-'0';
 bump9_V=[ones(size(bump9_bits,1),1),2*bump9_bits-1];
 bump9_pts=[zeros(1,bump9_n);.2,.1*ones(1,bump9_n-1);1,zeros(1,bump9_n-1);2*ones(1,bump9_n);3,zeros(1,bump9_n-1)];
 for bump9_i=1:size(bump9_pts,1)
  bump9_x=bump9_pts(bump9_i,:); bump9_sum=0;bump9_num=0;
  for bump9_j=1:size(bump9_V,1)
   bump9_v=bump9_V(bump9_j,:);bump9_d=dot(bump9_x,bump9_v);bump9_a=-bump9_d/bump9_n;bump9_A=1-dot(bump9_x,bump9_x)+bump9_d^2/bump9_n;
   if bump9_A<=0;continue;end
   bump9_h=sqrt(bump9_A/bump9_n);
   if bump9_a==0 || abs(abs(bump9_a)-bump9_h)<1e-13;bump9_log=0;else;bump9_log=bump9_a*(bump9_h^2-bump9_a^2)*log(abs((bump9_a+bump9_h)/(bump9_a-bump9_h)));end
   bump9_sum=bump9_sum+4*bump9_n^2/pi*(2*bump9_a^2*bump9_h-4*bump9_h^3/3+bump9_log);
   bump9_num=bump9_num+integral(@(bump9_s)bump9_integrand(bump9_s,bump9_x,bump9_v),0,Inf,'AbsTol',1e-9,'RelTol',1e-9)/pi;
  end
  assert(abs(bump9_sum-bump9_num)<1e-6);
  fprintf('BUMP9 n%d point%d formula %.17g second_diff %.17g error %.3g\n',bump9_n,bump9_i,bump9_sum,bump9_num,abs(bump9_sum-bump9_num));
 end
end
function bump9_val=bump9_integrand(bump9_s,bump9_x,bump9_v)
 bump9_q=dot(bump9_x,bump9_x);bump9_d=dot(bump9_x,bump9_v);bump9_n=dot(bump9_v,bump9_v);
 bump9_plus=1-bump9_q-2*bump9_d*bump9_s-bump9_n*bump9_s.^2;
 bump9_minus=1-bump9_q+2*bump9_d*bump9_s-bump9_n*bump9_s.^2;
 bump9_val=(max(bump9_plus,0).^2+max(bump9_minus,0).^2-2*max(1-bump9_q,0)^2)./bump9_s.^2;
 bump9_mask=(bump9_plus>0)&(bump9_minus>0);
 bump9_val(bump9_mask)=-4*bump9_n*(1-bump9_q)+8*bump9_d^2+2*bump9_n^2*bump9_s(bump9_mask).^2;
end
