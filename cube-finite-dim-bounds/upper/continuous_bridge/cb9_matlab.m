% cb9 fixed 18-point numerical cross-check; no search/clear/path changes.
cb9_root=fileparts(mfilename('fullpath'));
cb9_rec=jsondecode(fileread(fullfile(cb9_root,'cb9_results.json')));
cb9_cc=cell(20,1);cb9_count=0;
for cb9_i=1:3
    cb9_other=setdiff(1:3,cb9_i);cb9_j=cb9_other(1);cb9_k=cb9_other(2);
    for cb9_s=[-1,1]
        cb9_a=zeros(1,3);cb9_a(cb9_i)=cb9_s;cb9_a(cb9_j)=1;
        cb9_c=cb9_a;cb9_c(cb9_j)=-1;
        cb9_b=zeros(1,3);cb9_b(cb9_i)=cb9_s;cb9_b(cb9_k)=1;
        cb9_d=cb9_b;cb9_d(cb9_k)=-1;
        cb9_count=cb9_count+1;cb9_cc{cb9_count}=[cb9_a;cb9_b;cb9_c];
        cb9_count=cb9_count+1;cb9_cc{cb9_count}=[cb9_a;cb9_c;cb9_d];
    end
end
for cb9_s1=[-1,1]
 for cb9_s2=[-1,1]
  for cb9_s3=[-1,1]
   cb9_count=cb9_count+1;cb9_cc{cb9_count}=[cb9_s1 cb9_s2 0;cb9_s1 0 cb9_s3;0 cb9_s2 cb9_s3];
  end
 end
end
assert(cb9_count==20);
for cb9_i=1:20,assert(abs(det(cb9_cc{cb9_i}))==2);end
cb9_maxprel=0;cb9_maxgrel=0;
for cb9_idx=1:numel(cb9_rec.records)
    if iscell(cb9_rec.records),cb9_row=cb9_rec.records{cb9_idx};else,cb9_row=cb9_rec.records(cb9_idx);end
    cb9_x=zeros(1,3);
    for cb9_j=1:3
        cb9_parts=split(string(cb9_row.x_exact{cb9_j}),'/');cb9_x(cb9_j)=str2double(cb9_parts(1));
        if numel(cb9_parts)==2,cb9_x(cb9_j)=cb9_x(cb9_j)/str2double(cb9_parts(2));end
    end
    cb9_t=3/7;cb9_pcone=0;
    for cb9_i=1:20,cb9_pcone=cb9_pcone+1/prod(4*cb9_t-1i*(cb9_cc{cb9_i}*cb9_x.'));end
    cb9_pcone=real(cb9_pcone)/(4*pi^3);cb9_prat=cb9_heat(cb9_t,cb9_x);
    cb9_prel=abs(cb9_pcone-cb9_prat)/cb9_prat;assert(cb9_prel<1e-10);cb9_maxprel=max(cb9_maxprel,cb9_prel);
    cb9_g=cb9_green(cb9_x);
    if isinf(cb9_g)
        assert(strcmp(cb9_row.label,'origin') || strcmp(cb9_row.label,'body_diagonal'));
    else
        cb9_gref=str2double(string(cb9_row.G_mp60));cb9_grel=abs(cb9_g-cb9_gref)/cb9_gref;assert(cb9_grel<1e-9);cb9_maxgrel=max(cb9_maxgrel,cb9_grel);
        if strcmp(cb9_row.label,'axis'),assert(abs(cb9_g-3/(8*pi^3))<1e-12);end
        if strcmp(cb9_row.label,'coordinate_plane_edge'),assert(abs(cb9_g-(4*log(2)-1)/(8*pi^3))<1e-12);end
    end
end
fprintf('cb9 PASS: 20 det2 cones, 18 fixed points, max p relative diff %.17g, max G relative diff %.17g; axis/edge anchors and true infinite lines retained.\n',cb9_maxprel,cb9_maxgrel);
% Local helper functions are scoped to this portable script.
function cb9_p=cb9_heat(cb9_t,cb9_x)
 cb9_scale=max([4*cb9_t,abs(cb9_x)]);cb9_x=cb9_x/cb9_scale;cb9_a=4*cb9_t/cb9_scale;cb9_u=cb9_a^2;
 cb9_X=cb9_x(1)^2;cb9_Y=cb9_x(2)^2;cb9_Z=cb9_x(3)^2;cb9_S=cb9_X+cb9_Y+cb9_Z;cb9_S2=cb9_X*cb9_Y+cb9_Y*cb9_Z+cb9_Z*cb9_X;
 cb9_C=cb9_X*(cb9_Y-cb9_Z)^2+cb9_Y*(cb9_Z-cb9_X)^2+cb9_Z*(cb9_X-cb9_Y)^2;
 cb9_L=[cb9_x(1)-cb9_x(2),cb9_x(1)+cb9_x(2),cb9_x(1)-cb9_x(3),cb9_x(1)+cb9_x(3),cb9_x(2)-cb9_x(3),cb9_x(2)+cb9_x(3)];
 cb9_N=cb9_u*(20*cb9_u^3+24*cb9_S*cb9_u^2+4*(cb9_S^2+5*cb9_S2)*cb9_u+4*cb9_C);
 cb9_p=cb9_a*cb9_N/prod(cb9_u+cb9_L.^2)/(4*pi^3*cb9_scale^3);
end
function cb9_g=cb9_green(cb9_x)
 cb9_scale=max(abs(cb9_x));if cb9_scale==0,cb9_g=Inf;return;end
 cb9_x=sort(abs(cb9_x/cb9_scale),'descend');cb9_X=cb9_x(1)^2;cb9_Y=cb9_x(2)^2;cb9_Z=cb9_x(3)^2;
 cb9_S=cb9_X+cb9_Y+cb9_Z;cb9_T=cb9_S^2+5*(cb9_X*cb9_Y+cb9_Y*cb9_Z+cb9_Z*cb9_X);
 cb9_C=cb9_X*(cb9_Y-cb9_Z)^2+cb9_Y*(cb9_Z-cb9_X)^2+cb9_Z*(cb9_X-cb9_Y)^2;
 cb9_A=[(cb9_x(1)-cb9_x(2))^2,(cb9_x(1)+cb9_x(2))^2,(cb9_x(1)-cb9_x(3))^2,(cb9_x(1)+cb9_x(3))^2,(cb9_x(2)-cb9_x(3))^2,(cb9_x(2)+cb9_x(3))^2];
 cb9_zero=nnz(cb9_A==0);if cb9_zero>=3,cb9_g=Inf;return;end
 cb9_pts=unique(cb9_A(cb9_A>0)./(1+cb9_A(cb9_A>0)));
 cb9_g=integral(@(cb9_w)cb9_integrand(cb9_w,cb9_S,cb9_T,cb9_C,cb9_A,cb9_zero),0,1,'Waypoints',cb9_pts,'AbsTol',2e-11,'RelTol',2e-11)/(32*pi^3*cb9_scale^2);
end
function cb9_v=cb9_integrand(cb9_w,cb9_S,cb9_T,cb9_C,cb9_A,cb9_zero)
 cb9_z=1-cb9_w;cb9_den=ones(size(cb9_w));
 for cb9_a=cb9_A(cb9_A>0),cb9_den=cb9_den.*(cb9_w+cb9_a*cb9_z);end
 if cb9_zero==0
  cb9_num=cb9_w.*(20*cb9_w.^3+24*cb9_S*cb9_w.^2.*cb9_z+4*cb9_T*cb9_w.*cb9_z.^2+4*cb9_C*cb9_z.^3);
 elseif cb9_zero==1
  cb9_num=20*cb9_w.^3+24*cb9_S*cb9_w.^2.*cb9_z+4*cb9_T*cb9_w.*cb9_z.^2+4*cb9_C*cb9_z.^3;
 else
  cb9_num=20*cb9_w.^2+24*cb9_S*cb9_w.*cb9_z+4*cb9_T*cb9_z.^2;
 end
 cb9_v=cb9_num./cb9_den;
end
