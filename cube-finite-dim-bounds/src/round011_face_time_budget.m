% Actual scalar cross-check only, no proof claim.
ft11_n=[3,4];ft11_C=[256/3,32];ft11_A=[16/3,8];ft11_target=[.00057,.00019];
for ft11_j=1:2
 for ft11_eps=[.001,.0005]
  ft11_d=ft11_n(ft11_j);ft11_q=ft11_d/(ft11_d-1);
  ft11_E=2^(-(ft11_d-1)^2/ft11_d)*ft11_d/(ft11_d+1)*(ft11_C(ft11_j)/pi^ft11_d)^(1/ft11_d)*(ft11_A(ft11_j)/pi)^((ft11_d-1)/ft11_d)*ft11_eps^(1+1/ft11_d);
  ft11_B=(ft11_target(ft11_j)^(1/ft11_q)-ft11_E)^ft11_q;
  assert(ft11_B>0);fprintf('face_time n%d epsilon%.6g E=%.17g allowed_H=%.17g\n',ft11_d,ft11_eps,ft11_E,ft11_B);
 end
end
