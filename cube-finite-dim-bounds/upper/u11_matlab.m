% U11 four-dimensional body frame: minimal independent calibration.
% No cost search, no clear/path; preserve shared session.
upU11A=[1;0;0;0];upU11B=[1;1;0;0];upU11C=[1;1;1;0];
upU11D=[1;1;1;1];upU11E=[2;1;1;0];upU11F=[3;1;1;1];
upU11base=cat(3,[upU11A upU11B upU11E upU11F], ...
    [upU11E upU11B upU11C upU11D],[upU11E upU11B upU11D upU11F]);
upU11dets=[1 1 2];
for upU11k=1:3
    assert(abs(abs(det(upU11base(:,:,upU11k)))-upU11dets(upU11k))<1e-12);
end
upU11V=[];
for upU11i=0:7
    upU11V(end+1,:)=[1,2*(dec2bin(upU11i,3)-'0')-1];
end
upU11perms=perms(1:4);upU11x=[152/151,315/157,490/163,836/167];
upU11ss=0;upU11vol=0;upU11count=0;upU11p0sum=0;
for upU11i=1:24
 for upU11j=0:15
    upU11sg=2*(dec2bin(upU11j,4)-'0')-1;
    for upU11k=1:3
      upU11R=diag(upU11sg)*upU11base(upU11perms(upU11i,:),:,upU11k);
      upU11q=sum(abs(upU11V*upU11R),1);
      upU11b=(upU11x*upU11R)./upU11q;
      upU11alpha=upU11dets(upU11k)/prod(upU11q);
      upU11vol=upU11vol+upU11alpha/24;
      upU11p0sum=upU11p0sum+upU11alpha;
      if all(upU11b>0) || all(upU11b<0)
        upU11dd=0;
      else
        upU11dd=0;
        for upU11h=1:4
          upU11other=setdiff(1:4,upU11h);
          assert(all(abs(upU11b(upU11h)-upU11b(upU11other))>1e-12));
          upU11dd=upU11dd+sign(upU11b(upU11h))/prod(upU11b(upU11h)-upU11b(upU11other));
        end
      end
      upU11ss=upU11ss+upU11alpha*upU11dd;
      upU11count=upU11count+1;
    end
 end
end
upU11G=-upU11ss/(32*pi^3);
upU11Gexact=0.0000138035907751128585307687051;
assert(upU11count==1152 && abs(upU11vol-11/4608)<1e-14);
assert(abs(upU11G-upU11Gexact)/upU11Gexact<1e-8);
upU11I=integral(@(upU11t)1./((upU11t.^2+1).*(upU11t.^2+4)),0,Inf);
upU11Irep=integral(@(upU11t)1./(upU11t.^2+1).^2,0,Inf);
upU11eps=1/8;
upU11Izero=integral(@(upU11t)1./(upU11t.^2.*(upU11t.^2+1)),upU11eps,Inf);
assert(abs(upU11I-pi/12)<1e-12 && abs(upU11Irep-pi/4)<1e-12);
assert(abs(upU11Izero-(1/upU11eps-pi/2+atan(upU11eps)))<1e-12);
fprintf('U11 cone_count %d spectral_volume %.17g expected %.17g\n',upU11count,upU11vol,11/4608);
fprintf('U11 p_t0_pi4_coefficient %.17g expected %.17g\n',upU11p0sum/16,11/3072);
fprintf('U11 time_pi12 %.17g defect %.17g repeated_pi4_defect %.17g zero_cutoff_defect %.17g\n',upU11I,abs(upU11I-pi/12),abs(upU11Irep-pi/4),abs(upU11Izero-(1/upU11eps-pi/2+atan(upU11eps))));
fprintf('U11 generic_first_G %.17g exact_display %.17g relative_defect %.17g\n',upU11G,upU11Gexact,abs(upU11G-upU11Gexact)/upU11Gexact);
fprintf('U11 scope: representation calibration only; no all-face minimum or cost.\n');
