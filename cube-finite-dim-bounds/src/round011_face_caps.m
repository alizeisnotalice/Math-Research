% Main serial numerical crosscheck of exact finite determinant certificate.
fc11_ns=[3 3 3 3 4];fc11_taus=[.1 1/3 1 3 0];
fc11_expected=[5892 36800/3 38400 187200 5632];
for fc11_case=1:5
 fc11_n=fc11_ns(fc11_case);fc11_tau=fc11_taus(fc11_case);
 fc11_body=[ones(2^(fc11_n-1),1),2*(dec2bin(0:2^(fc11_n-1)-1)-'0')-1];
 fc11_V=fc11_body;fc11_w=ones(size(fc11_body,1),1);
 if fc11_tau>0
  for fc11_i=1:fc11_n
   for fc11_j=fc11_i+1:fc11_n
    for fc11_sign=[-1 1]
     fc11_v=zeros(1,fc11_n);fc11_v(fc11_i)=1;fc11_v(fc11_j)=fc11_sign;
     fc11_V=[fc11_V;fc11_v];fc11_w=[fc11_w;fc11_tau];
    end
   end
  end
 end
 fc11_m=size(fc11_V,1);fc11_combos=nchoosek(1:size(fc11_body,1),fc11_n-2);fc11_S=0;
 for fc11_i=1:fc11_m
  for fc11_j=1:fc11_m
   if fc11_i==fc11_j,continue,end
   fc11_det=0;
   for fc11_b=1:size(fc11_combos,1)
    fc11_d=abs(det([fc11_V(fc11_i,:);fc11_V(fc11_j,:);fc11_body(fc11_combos(fc11_b,:),:)]));
    assert(abs(fc11_d-round(fc11_d))<1e-10);fc11_det=max(fc11_det,round(fc11_d));
   end
   assert(fc11_det>0);fc11_S=fc11_S+fc11_w(fc11_i)*fc11_w(fc11_j)/fc11_det;
  end
 end
 fc11_Cnum=4*fc11_m^2*fc11_n*fc11_S;
 assert(abs(fc11_Cnum-fc11_expected(fc11_case))<1e-7);
 fprintf('11A PASS n=%d tau=%.12g pairs=%d C*pi^n=%.14g\n',fc11_n,fc11_tau,fc11_m*(fc11_m-1),fc11_Cnum);
end
