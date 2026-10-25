% cb8 signed certificate numeric cross-check; no solver/clear/path changes.
cb8_root=fileparts(mfilename('fullpath'));
cb8_cases=[3 1;3 2;4 1];
cb8_BN=[0 0 0 1;1 2 4 4;4 4 2 1;1 0 0 0];
cb8_DN=[0 1 2 3;3 0 -3 -6;-6 -3 0 3;3 2 1 0];
for cb8_case=1:3
    cb8_n=cb8_cases(cb8_case,1);cb8_q=cb8_cases(cb8_case,2);
    cb8_rec=jsondecode(fileread(fullfile(cb8_root,sprintf('cb8_certificate_n%d_q%d.json',cb8_n,cb8_q))));
    cb8_K=cb8_rec.K;cb8_den=str2double(cb8_rec.common_coefficient_denominator);
    cb8_a=str2double(cb8_rec.common_coefficient_numerators)/cb8_den;
    cb8_points=cb8_rec.orbits;cb8_shape=repmat(4,1,cb8_n);
    cb8_origin=sum(cb8_rec.origin_integer_row(:).*cb8_a(:))/3^cb8_n;
    cb8_hdiag=0;
    for cb8_j=1:size(cb8_points,1)
        if max(cb8_points(cb8_j,:))<=1
            cb8_s=sum(cb8_points(cb8_j,:)==1);
            cb8_hdiag=cb8_hdiag+cb8_q^2*cb8_rec.origin_integer_row(cb8_j)*cb8_a(cb8_j)*(9*cb8_s/cb8_n-3)/3^cb8_n;
        end
    end
    cb8_grids=cell(1,cb8_n);[cb8_grids{:}]=ndgrid(0:cb8_K+1);
    cb8_all=zeros(numel(cb8_grids{1}),cb8_n);
    for cb8_j=1:cb8_n,cb8_all(:,cb8_j)=cb8_grids{cb8_j}(:);end
    cb8_cells=cb8_all(all(diff(cb8_all,1,2)>=0,2),:);
    cb8_minpsi=inf;cb8_maxlap=-inf;cb8_cubeint=0;cb8_fullcount=0;
    for cb8_c=1:size(cb8_cells,1)
        cb8_cell=cb8_cells(cb8_c,:);cb8_choices=cell(1,cb8_n);
        for cb8_j=1:cb8_n,cb8_choices{cb8_j}=max(-cb8_K,cb8_cell(cb8_j)-1):min(cb8_K,cb8_cell(cb8_j)+2);end
        [cb8_grids{:}]=ndgrid(cb8_choices{:});cb8_centers=zeros(numel(cb8_grids{1}),cb8_n);
        for cb8_j=1:cb8_n,cb8_centers(:,cb8_j)=cb8_grids{cb8_j}(:);end
        [cb8_tf,cb8_idx]=ismember(sort(abs(cb8_centers),2),cb8_points,'rows');assert(all(cb8_tf));
        cb8_poly=zeros(cb8_shape);cb8_lap=zeros(cb8_shape);
        for cb8_k=1:size(cb8_centers,1)
            cb8_coeff=cb8_a(cb8_idx(cb8_k));if cb8_coeff==0,continue;end
            cb8_m=cb8_cell-cb8_centers(cb8_k,:)+3;
            cb8_b=cell(1,cb8_n);cb8_d=cell(1,cb8_n);
            for cb8_j=1:cb8_n
                cb8_sh=ones(1,cb8_n);cb8_sh(cb8_j)=4;
                cb8_b{cb8_j}=reshape(cb8_BN(cb8_m(cb8_j),:),cb8_sh);
                cb8_d{cb8_j}=reshape(cb8_DN(cb8_m(cb8_j),:),cb8_sh);
            end
            cb8_term=ones(cb8_shape);
            for cb8_j=1:cb8_n,cb8_term=cb8_term.*cb8_b{cb8_j};end
            cb8_poly=cb8_poly+cb8_coeff*cb8_term/6^cb8_n;
            for cb8_ax=1:cb8_n
                cb8_term=ones(cb8_shape);
                for cb8_j=1:cb8_n
                    if cb8_j==cb8_ax,cb8_term=cb8_term.*cb8_d{cb8_j};else,cb8_term=cb8_term.*cb8_b{cb8_j};end
                end
                cb8_lap=cb8_lap+cb8_coeff*cb8_q^2*cb8_term/(3*6^(cb8_n-1));
            end
        end
        cb8_mult=2^cb8_n*size(unique(perms(cb8_cell),'rows'),1);
        cb8_fullcount=cb8_fullcount+cb8_mult;
        cb8_minpsi=min(cb8_minpsi,min(cb8_poly(:)));cb8_maxlap=max(cb8_maxlap,max(cb8_lap(:)));
        if max(cb8_cell)<cb8_q,cb8_cubeint=cb8_cubeint+cb8_mult*sum(cb8_lap(:))/(4^cb8_n*cb8_q^cb8_n);end
    end
    cb8_value=1-cb8_cubeint/2^cb8_n;
    cb8_parts=split(string(cb8_rec.D_exact),'/');cb8_expected=str2double(cb8_parts(1));if numel(cb8_parts)==2,cb8_expected=cb8_expected/str2double(cb8_parts(2));end
    assert(cb8_fullcount==(2*cb8_K+4)^cb8_n && cb8_minpsi>=-1e-8 && cb8_maxlap<=1+1e-8);
    assert(abs(cb8_origin)<1e-8 && abs(cb8_value-cb8_expected)<1e-8 && cb8_hdiag>=-1e-8 && cb8_n*cb8_hdiag<=1+1e-8);
    fprintf('cb8 n=%d q=%d PASS: %d covered full cells (%d evaluated representatives), minPsi %.17g, maxDelta %.17g, psi0 %.17g, Hdiag %.17g, D %.17g\n',cb8_n,cb8_q,cb8_fullcount,size(cb8_cells,1),cb8_minpsi,cb8_maxlap,cb8_origin,cb8_hdiag,cb8_value);
end
% Floating cross-check only. Exact positivity/origin are from rational certificates.
