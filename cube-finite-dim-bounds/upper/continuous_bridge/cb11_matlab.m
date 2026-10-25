% cb11 36 preregistered generic parameter cases; no cost/min search.
cb11_root=fileparts(mfilename('fullpath'));
cb11_data=jsondecode(fileread(fullfile(cb11_root,'cb11_results.json')));
cb11_bases={ [1 0 0;1 1 0;2 1 1], [1 1 0;2 1 1;1 1 1] };
cb11_pp=perms(1:3);cb11_ss=2*(dec2bin(0:7,3)-'0')-1;cb11_cones=cell(96,1);cb11_count=0;
for cb11_ip=1:6
 for cb11_is=1:8
  for cb11_ib=1:2
   cb11_count=cb11_count+1;cb11_cones{cb11_count}=cb11_bases{cb11_ib}(:,cb11_pp(cb11_ip,:)).*cb11_ss(cb11_is,:);assert(abs(det(cb11_cones{cb11_count}))==1);
  end
 end
end
cb11_maxrel=0;cb11_maxprel=0;
for cb11_idx=1:numel(cb11_data.records)
 if iscell(cb11_data.records),cb11_row=cb11_data.records{cb11_idx};else,cb11_row=cb11_data.records(cb11_idx);end
 cb11_part=split(string(cb11_row.tau_exact),'/');cb11_tau=str2double(cb11_part(1));if numel(cb11_part)==2,cb11_tau=cb11_tau/str2double(cb11_part(2));end
 cb11_x=zeros(1,3);
 for cb11_j=1:3
  cb11_part=split(string(cb11_row.x_exact{cb11_j}),'/');cb11_x(cb11_j)=str2double(cb11_part(1));if numel(cb11_part)==2,cb11_x(cb11_j)=cb11_x(cb11_j)/str2double(cb11_part(2));end
 end
 cb11_G=0;cb11_p=0;cb11_peak=0;
 for cb11_i=1:96
  cb11_R=cb11_cones{cb11_i};cb11_sorted=sort(abs(cb11_R),2,'descend');
  cb11_rates=max(4*cb11_sorted(:,1),2*sum(cb11_sorted,2))+2*cb11_tau*(2*cb11_sorted(:,1)+cb11_sorted(:,2));
  cb11_b=(cb11_R*cb11_x.')./cb11_rates;cb11_dd=0;
  for cb11_j=1:3
   cb11_other=setdiff(1:3,cb11_j);cb11_dd=cb11_dd+log(abs(cb11_b(cb11_j)))/prod(cb11_b(cb11_j)-cb11_b(cb11_other));
  end
  cb11_G=cb11_G+cb11_dd/prod(cb11_rates);
  cb11_p=cb11_p+1/prod((3/7)*cb11_rates-1i*(cb11_R*cb11_x.'));
  cb11_peak=cb11_peak+1/prod(cb11_rates);
 end
 cb11_G=cb11_G/(8*pi^3);cb11_p=real(cb11_p)/(8*pi^3);cb11_peak=cb11_peak/8;
 cb11_ref=str2double(string(cb11_row.G_mp80));cb11_rel=abs(cb11_G-cb11_ref)/cb11_ref;cb11_prel=abs(cb11_p-cb11_row.p_full96_complex)/cb11_row.p_full96_complex;
 assert(cb11_rel<1e-8 && cb11_prel<1e-10);
 assert(abs(cb11_peak-5/(8*(1+cb11_tau)*(2+3*cb11_tau)*(4+5*cb11_tau)))<1e-12);
 cb11_maxrel=max(cb11_maxrel,cb11_rel);cb11_maxprel=max(cb11_maxprel,cb11_prel);
end
fprintf('cb11 PASS: 96 det1 cones, 36 fixed parameter points, max Green relative diff %.17g, max heat relative diff %.17g; peak coefficients verified.\n',cb11_maxrel,cb11_maxprel);
% Numeric generic calibration only. Zero/repeated knots use exact rules in Python/proof.
