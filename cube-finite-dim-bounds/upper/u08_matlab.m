% U08: the same one matrix and one analytic family; main runs serially.
upU08N=5;upU08t=1/4;upU08Q=zeros(5);
for upU08i=1:5
    for upU08j=1:5
        if upU08i~=upU08j
            upU08d=min(abs(upU08i-upU08j),5-abs(upU08i-upU08j));
            upU08Q(upU08i,upU08j)=1/upU08d^2;
        end
    end
    upU08Q(upU08i,upU08i)=-sum(upU08Q(upU08i,:));
end
upU08P=eye(5)+upU08t*upU08Q;upU08f=[2;0;0;0;0];upU08u=[2/5;0;0;0;0];
upU08sigma=upU08f+(upU08P*upU08u-upU08u)/upU08t;
upU08energy=upU08u'*(upU08u-upU08P*upU08u)/upU08t;
fprintf('U08 P row defect %.3g fixedpoint defect %.3g mass %.17g energy %.17g\n', ...
    max(abs(sum(upU08P,2)-1)),max(abs(max(0,upU08P*upU08u+upU08t*(upU08f-1))-upU08u)),sum(upU08sigma),upU08energy);
fprintf('U08 sigma %.17g %.17g %.17g %.17g %.17g\n',upU08sigma);
upU08w=zeros(5,1);
for upU08j=1:20
    upU08w=max(0,upU08P*upU08w+upU08t*([6;0;0;0;0]-1));
end
fprintf('U08 supercritical 20-step mass %.17g (proved lower bound 5)\n',sum(upU08w));
for upU08n=[3,4]
    upU08sourcePeak=100/pi^upU08n;upU08sigmaPeak=100/(4*pi)^upU08n;
    upU08uPeak=100/(pi^upU08n*(upU08n-1))*(1-4^(1-upU08n));
    fprintf('U08 analytic n%d f_peak %.17g sigma_peak %.17g uT_peak %.17g uT_L1 300; saturation fails\n', ...
        upU08n,upU08sourcePeak,upU08sigmaPeak,upU08uPeak);
end
upU08conv=integral(@(upU08y) (.5./(pi*(.5^2+(.7-upU08y).^2))).*(1./(pi*(1+upU08y.^2))),-Inf,Inf,'AbsTol',1e-12,'RelTol',1e-11);
upU08exact=1.5/(pi*(1.5^2+.7^2));
fprintf('U08 Cauchy semigroup one-dimensional convolution %.17g exact %.17g defect %.3g\n',upU08conv,upU08exact,abs(upU08conv-upU08exact));
