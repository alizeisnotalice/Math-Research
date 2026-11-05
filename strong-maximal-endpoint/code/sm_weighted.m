function rows=sm_weighted
% Hypoexponential cube-volume controls; fixed dimension then depth->infinity.
root=fileparts(fileparts(mfilename('fullpath')));
rows=cell(0,1);
for weights={[1 0.5 0.25],[1 1 0.5]}
    a=weights{1}; rates=1./a; n=numel(a);
    Q=diag(-rates)+diag(rates(1:end-1),1);
    m=sum(a==1); coefficient=prod(1./(1-a(a<1)))/factorial(m-1);
    for L=[2 4 8 16 32]
        tail=[1 zeros(1,n-1)]*expm(Q*L)*ones(n,1);
        density=[1 zeros(1,n-1)]*expm(Q*L)*(-Q)*ones(n,1);
        scaled=exp(L)*tail/L^(m-1);
        row=struct('weights',a,'L',L,'tail',tail,'density',density, ...
            'scaled_tail',scaled,'asymptotic_coefficient',coefficient);
        rows{end+1}=row; %#ok<AGROW>
        fprintf('weights=%s L=%g scaledTail=%.12g target=%.12g\n', ...
            mat2str(a),L,scaled,coefficient);
    end
end
fid=fopen(fullfile(root,'rounds','r006','weighted.json'),'w');
assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(struct('matlab_version',version,'rows',{rows}),PrettyPrint=true));
end
