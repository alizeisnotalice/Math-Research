% Independent strict-feasibility LP, then 2D clipping / 3D vertices.
% Floating point volume checks are not proof certificates.
base=fileparts(mfilename('fullpath'));started=tic;
parameters=[2 .5;2 1;2 1.5;3 .5;3 1.5;3 2.5];
amplitudes=[0 .1 .25 .5];rows=zeros(24,6);count=0;
options=optimoptions('linprog','Display','none','Algorithm','dual-simplex');
for j=1:size(parameters,1)
    n=parameters(j,1);beta=parameters(j,2);
    for amplitude=amplitudes
        h=amplitude*linspace(-1,1,n);f=reshape([h;h],[],1);
        masks=zeros(2^(2*n),2*n);
        for bit=1:2*n,masks(:,bit)=bitget(uint32((0:2^(2*n)-1)'),bit);end
        u=zeros(size(masks,1),n);
        for a=1:n,u(:,a)=-(masks(:,2*a-1)-masks(:,2*a));end
        q=n-beta-sum(masks,2)/2;
        rhs=1+(1-beta/n)*sum(f)-masks*f;
        % Strict feasibility is checked BEFORE discarding zero projected normals.
        [solution,~,flag]=linprog([zeros(n+1,1);-1],[u q ones(size(q))],rhs, ...
            [],[],[],[],options);
        assert(flag>0);
        margin=solution(end);
        if margin<=1e-9
            volume=0;constraintError=0;
        else
            [a,b]=projectScale(u,q,rhs);
            if n==2
                radius=(1+2*sum(abs(f)))/min(1,2*min(beta,n-beta))+1;
                polygon=radius*[-1 -1;1 -1;1 1;-1 1];
                for k=1:size(a,1),polygon=clipPolygon(polygon,a(k,:),b(k));end
                volume=polyarea(polygon(:,1),polygon(:,2));vertices=polygon;
            else
                % Remove redundant facets independently by linear programming.
                k=1;
                while k<=size(a,1)
                    keep=true(size(b));keep(k)=false;
                    [~,value,status]=linprog(-a(k,:)',a(keep,:),b(keep),[],[],[],[],options);
                    if status>0 && -value<=b(k)+2e-10
                        a(k,:)=[];b(k)=[];
                    else
                        k=k+1;
                    end
                end
                combinations=nchoosek(1:size(a,1),3);
                vertices=zeros(size(combinations,1),3);nv=0;
                for k=1:size(combinations,1)
                    matrix=a(combinations(k,:),:);
                    if rcond(matrix)<1e-11,continue;end
                    point=matrix\b(combinations(k,:));
                    if all(a*point<=b+2e-9),nv=nv+1;vertices(nv,:)=point';end
                end
                vertices=unique(round(vertices(1:nv,:),12),'rows');
                assert(size(vertices,1)>=4);[~,volume]=convhulln(vertices);
            end
            constraintError=max(a*vertices'-b,[],'all');
            assert(constraintError<3e-8);
        end
        count=count+1;rows(count,:)=[n beta amplitude volume margin constraintError];
    end
end
assert(rows(rows(:,1)==2 & rows(:,2)==1 & rows(:,3)==.5,4)==0);
result=struct('kind','independent-strict-LP-and-geometry-numerical-evidence', ...
    'version',version,'precision','IEEE754 double','elapsed_seconds',toc(started), ...
    'strict_margin_guard',1e-9,'facet_redundancy_tolerance',2e-10, ...
    'vertex_feasibility_tolerance',2e-9,'random_seed',[], ...
    'columns',{{'n','beta','amplitude','volume','strict_extended_margin','constraint_residual'}}, ...
    'rows',rows);
fid=fopen(fullfile(base,'..','results','global_leading_matlab.json'),'w');
assert(fid~=-1);fprintf(fid,'%s\n',jsonencode(result));fclose(fid);
fprintf('Strict face fibres: %d volumes in %.3f seconds; flat threshold excluded.\n',count,result.elapsed_seconds);

function [a,b]=projectScale(u,q,rhs)
    positive=find(q>1e-12);negative=find(q< -1e-12);zero=find(abs(q)<=1e-12);
    a=u(zero,:);b=rhs(zero);
    for p=positive'
        for m=negative'
            a(end+1,:)=q(p)*u(m,:)-q(m)*u(p,:); %#ok<AGROW>
            b(end+1,1)=q(p)*rhs(m)-q(m)*rhs(p); %#ok<AGROW>
        end
    end
    norms=vecnorm(a,2,2);assert(all(b(norms<1e-12)>0));
    a=a(norms>=1e-12,:);b=b(norms>=1e-12);norms=norms(norms>=1e-12);
    a=a./norms;b=b./norms;
    [~,first,group]=unique(round(a,12),'rows');
    a=a(first,:);b=accumarray(group,b,[],@min);
end

function result=clipPolygon(polygon,normal,offset)
    result=zeros(2*size(polygon,1),2);count=0;
    for j=1:size(polygon,1)
        p=polygon(j,:);q=polygon(mod(j,size(polygon,1))+1,:);
        dp=dot(normal,p)-offset;dq=dot(normal,q)-offset;
        if dp<=1e-12,count=count+1;result(count,:)=p;end
        if (dp<0 && dq>0)||(dp>0 && dq<0)
            count=count+1;result(count,:)=p+dp/(dp-dq)*(q-p);
        end
    end
    result=result(1:count,:);
end
