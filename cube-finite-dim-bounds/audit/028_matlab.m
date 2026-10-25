% Parent serial crosscheck only. Set aud028_project before eval(fileread(...)).
% No clear, path changes, LP, or Green routine. Numerical, not a proof.
assert(exist('aud028_project','var')==1,'Set aud028_project to the project directory.');
aud028_z=jsondecode(fileread(fullfile(aud028_project,'results','round010c_source_input.json')));
aud028_x=cellfun(@aud028_fraction,aud028_z.point).';
aud028_w=cellfun(@aud028_fraction,aud028_z.weights);
aud028_V=[1 -1 -1;1 -1 1;1 1 -1;1 1 1];
aud028_e=sum(aud028_x.^2);
aud028_c=(4*log(2)-1)/(24*pi^3);
aud028_R=15/(2*(4*log(2)-1));
assert(aud028_e<aud028_R);
aud028_mu=0;
for aud028_j=1:4
    aud028_d=aud028_x*aud028_V(aud028_j,:).';
    aud028_q=sqrt(aud028_d^2+3*(aud028_R-aud028_e));
    for aud028_sigma=[-1,1]
        aud028_U=(aud028_q+aud028_sigma*aud028_d)/(aud028_R-aud028_e);
        aud028_mu=aud028_mu+aud028_c/pi*integral(@(aud028_u) ...
            1-aud028_R*aud028_u.^2./(aud028_e*aud028_u.^2+ ...
            2*aud028_sigma*aud028_d*aud028_u+3),0,aud028_U, ...
            'AbsTol',1e-13,'RelTol',1e-12);
    end
end
aud028_contributions=zeros(8,1);aud028_counts=zeros(8,1);
aud028_signs=2*(dec2bin(0:7)-'0')-1;
for aud028_j=1:8
    aud028_seed=cellfun(@aud028_fraction,aud028_z.specs(aud028_j).seed).';
    aud028_radius=aud028_fraction(aud028_z.specs(aud028_j).radius);
    aud028_permutations=perms(aud028_seed);aud028_centers=[];
    for aud028_k=1:size(aud028_permutations,1)
        aud028_centers=[aud028_centers;aud028_signs.*aud028_permutations(aud028_k,:)]; %#ok<AGROW>
    end
    aud028_centers=unique(aud028_centers,'rows');
    aud028_counts(aud028_j)=size(aud028_centers,1);
    for aud028_k=1:size(aud028_centers,1)
        aud028_y=(aud028_x-aud028_centers(aud028_k,:))/aud028_radius;
        for aud028_l=1:4
            aud028_contributions(aud028_j)=aud028_contributions(aud028_j)+ ...
                aud028_w(aud028_j)/aud028_radius*aud028_second_difference(aud028_y,aud028_V(aud028_l,:));
        end
    end
end
aud028_LB=sum(aud028_contributions);
aud028_negative=aud028_mu-aud028_LB;
assert(isequal(aud028_counts,[1;6;1;8;24;12;6;1]));
assert(abs(aud028_mu-.00376150414280235821)<1e-11);
assert(abs(aud028_LB-.00578292705177662161)<1e-10);
assert(aud028_negative<-.00202142);
fprintf('028 PASS centers=%d mu_upper=%.17g LB=%.17g source_upper=%.17g\n', ...
    sum(aud028_counts),aud028_mu,aud028_LB,aud028_negative);
disp(aud028_contributions.');

function aud028_f=aud028_fraction(aud028_text)
    aud028_parts=split(string(aud028_text),'/');
    aud028_f=str2double(aud028_parts(1));
    if numel(aud028_parts)==2,aud028_f=aud028_f/str2double(aud028_parts(2));end
end

function aud028_value=aud028_second_difference(aud028_y,aud028_v)
    % Integral on s>0 of [b(y+sv)+b(y-sv)-2b(y)]/s^2.
    aud028_d=dot(aud028_y,aud028_v);aud028_k=1-dot(aud028_y,aud028_y);
    aud028_discriminant=aud028_d^2+3*aud028_k;
    if aud028_discriminant<=0,aud028_value=0;return;end
    aud028_roots=(-aud028_d+[-1,1]*sqrt(aud028_discriminant))/3;
    aud028_cuts=unique([0,abs(aud028_roots)]);
    aud028_f0=max(aud028_k,0)^2;aud028_value=0;
    for aud028_i=1:numel(aud028_cuts)-1
        aud028_a=aud028_cuts(aud028_i);aud028_b=aud028_cuts(aud028_i+1);
        aud028_mid=(aud028_a+aud028_b)/2;
        if aud028_k-2*aud028_d*aud028_mid-3*aud028_mid^2>0 && ...
                aud028_k+2*aud028_d*aud028_mid-3*aud028_mid^2>0
            % Stable exact polynomial on the interval containing s=0.
            aud028_value=aud028_value+(8*aud028_d^2-12*aud028_k)*(aud028_b-aud028_a) ...
                +6*(aud028_b^3-aud028_a^3);
        elseif aud028_b>aud028_a
            aud028_value=aud028_value+integral(@(aud028_s) ...
                (max(aud028_k-2*aud028_d*aud028_s-3*aud028_s.^2,0).^2+ ...
                 max(aud028_k+2*aud028_d*aud028_s-3*aud028_s.^2,0).^2-2*aud028_f0)./aud028_s.^2, ...
                aud028_a,aud028_b,'AbsTol',1e-11,'RelTol',1e-10);
        end
    end
    aud028_value=(aud028_value-2*aud028_f0/max(aud028_cuts))/pi;
end
