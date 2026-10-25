function cube5_round002_matlab()
tic; base=zeros(5,5,7);
base(:,:,1)=[1 1 4 3 2;1 0 1 1 1;0 0 1 1 1;0 0 1 1 0;0 0 1 0 0];
base(:,:,2)=[1 4 2 3 2;1 1 1 1 1;0 1 1 1 1;0 1 1 1 0;0 1 1 0 0];
base(:,:,3)=[1 1 2 3 2;1 1 1 1 1;1 0 1 1 1;1 0 1 1 0;0 0 1 0 0];
base(:,:,4)=[1 2 1 2 2;1 2 1 1 1;1 2 0 1 1;1 1 0 1 0;0 1 0 1 0];
base(:,:,5)=[1 2 1 3 2;1 2 1 3 1;1 2 0 2 1;1 1 0 2 1;0 1 0 2 1];
base(:,:,6)=[1 1 2 3 2;1 1 2 3 1;1 1 2 2 1;1 1 1 2 1;1 0 1 2 1];
base(:,:,7)=[1 1 2 1 2;1 1 2 1 1;1 1 2 0 1;1 0 1 0 0;0 0 1 0 0];
V=[ones(16,1),2*(dec2bin(0:15,4)-'0')-1];
P=perms(1:5); S=2*(dec2bin(0:31,5)-'0')-1;
x=[1,1/7,2/7,3/7,4/7]; total=0; counter=0;
for k=1:7
 R0=base(:,:,k); d=round(abs(det(R0))); rate=sum(abs(V*R0),1);
 for ip=1:size(P,1)
  Rperm=R0(P(ip,:),:);
  for is=1:32
   R=(S(is,:).'*ones(1,5)).*Rperm;
   total=total+d/prod(rate-1i*(x*R)); counter=counter+1;
  end
 end
end
value=total/(2*pi)^5;
fprintf('correct-column Fourier: terms=%d p1=%.17g imag=%.3g elapsed=%.6g s\n',counter,real(value),imag(value),toc);
fprintf('difference from Arb midpoint: %.3g\n',real(value)-1.47717170324864630922e-7);
end
