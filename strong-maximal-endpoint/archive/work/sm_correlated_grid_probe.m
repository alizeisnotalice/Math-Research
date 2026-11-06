function result = sm_correlated_grid_probe()
% Finite 2-D strong-rectangle probe.  All reported candidates are rechecked
% from the rational amplitude grid; MILP statuses are retained verbatim.
root = fileparts(mfilename('fullpath'));
jsonPath = fullfile(root,'sm_correlated_grid_probe.json');
notePath = fullfile(root,'sm_correlated_grid_probe_notes.md');
amp = [0 1 3/2 2 4 8 16]; phi = amp.*(1+max(log(max(amp,1)),0));
delta = 1/1000; inEdges = 0:4;
outEdges = [-16 -4 -1 0 1 2 3 4 5 8 20];
[rect, overlap, area, recvCover] = buildGrid(inEdges,outEdges);
tinyHandCheck(rect,overlap,area,recvCover,outEdges);
keep = max(amp)*sum(overlap,2) > (1+delta)*area;
rect = rect(keep,:); overlap = overlap(keep,:); area = area(keep); recvCover = recvCover(:,keep);
R=size(rect,1); C=16; L=numel(amp); W=size(recvCover,1);
% y(cell,level), z(rectangle), w(receiving cell)
ny=C*L; iz=ny+(1:R); iw=ny+R+(1:W); N=ny+R+W;
A=[]; b=[];
% Each selected rectangle has actual average at least 1+delta.
for q=1:R
    row=sparse(1,N); row(1:C*L)=reshape((-overlap(q,:)'*amp).',C*L,1)'; row(iz(q))=(1+delta)*area(q);
    A=[A;row]; b=[b;0]; %#ok<AGROW>
end
% A receiving cell can be counted only if a selected rectangle contains it.
for j=1:W
    row=sparse(1,N); row(iw(j))=1; row(iz)=-recvCover(j,:);
    A=[A;row]; b=[b;0]; %#ok<AGROW>
end
A=[A; -sparse(1,iw,ones(1,W),1,N)]; b=[b;-1];
Aeq=sparse(C,N); beq=ones(C,1);
for c=1:C, Aeq(c,(c-1)*L+(1:L))=1; end
lb=zeros(N,1); ub=ones(N,1); intcon=1:N;
options=optimoptions('intlinprog','Display','off','MaxTime',5,'RelativeGapTolerance',1e-6);
lambda=0; history=struct([]); best=[];
for iter=1:6
    f=zeros(N,1); f(1:ny)=repmat((lambda*phi)',C,1); f(iw)=-receivingAreas(outEdges)';
    [x,fval,exitflag,out]=intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,[],options);
    rec=struct('iteration',iter,'lambda_in',lambda,'exitflag',exitflag,'fval',fval, ...
        'message',fieldOr(out,'message',''),'relativegap',fieldOr(out,'relativegap',NaN), ...
        'numnodes',fieldOr(out,'numnodes',NaN),'constrviolation',fieldOr(out,'constrviolation',NaN));
    if isempty(x), history(end+1)=rec; break; end %#ok<AGROW>
    candidate = recheck(round(x(1:ny)),round(x(iz)),amp,phi,overlap,area,recvCover,outEdges,delta,rect);
    rec.halo_area=candidate.halo_area; rec.cost=candidate.cost; rec.ratio=candidate.ratio;
    rec.dinkelbach_residual=candidate.halo_area-lambda*candidate.cost;
    rec.recheck_ok=candidate.recheck_ok;
    if isempty(history), history=rec; else, history(end+1)=rec; end %#ok<AGROW>
    if isempty(best) || candidate.ratio>best.ratio, best=candidate; end
    fprintf('iteration=%d halo=%g cost=%g ratio=%g exit=%d\n',iter,candidate.halo_area,candidate.cost,candidate.ratio,exitflag);
    lambda=candidate.ratio;
    if exitflag>0 && abs(rec.dinkelbach_residual)<1e-8, break; end
end
baseline=axisBlockBaseline(amp,phi,overlap,area,recvCover,outEdges,delta);
result=struct('model','finite 4x4 unit input cells; real axis-parallel grid-edge rectangles; binary amplitude one-hot, rectangle, and receiving-cell variables', ...
    'amplitudes',amp,'orlicz_costs',phi,'delta',delta,'out_edges',outEdges, ...
    'rectangle_count_after_prune',R,'solver_history',history,'best_candidate',best, ...
    'axisblock_constant_baseline',baseline,'interpretation','finite-grid family probe only; a time-limit incumbent is not recorded as an optimizer');
writeText(jsonPath,jsonencode(result,PrettyPrint=true));
writeText(notePath,notes(result));
fprintf('rectangles=%d best ratio=%.12g baseline=%.12g\n',R,best.ratio,baseline.ratio);
end

function [rect,overlap,area,cover]=buildGrid(inE,outE)
R=[]; O=[]; A=[]; C=[]; [gx,gy]=ndgrid(1:numel(outE)-1,1:numel(outE)-1);
for ax=1:numel(outE)-1, for bx=ax+1:numel(outE)
for ay=1:numel(outE)-1, for by=ay+1:numel(outE)
    rr=[outE(ax) outE(bx) outE(ay) outE(by)]; R=[R;rr]; %#ok<AGROW>
    A=[A;(rr(2)-rr(1))*(rr(4)-rr(3))]; %#ok<AGROW>
    oo=zeros(1,16); k=0;
    for i=1:4, for j=1:4
        k=k+1; oo(k)=max(0,min(i,rr(2))-max(i-1,rr(1)))*max(0,min(j,rr(4))-max(j-1,rr(3)));
    end,end
    O=[O;oo]; %#ok<AGROW>
    inside=(outE(gx(:))>=rr(1)) & (outE(gx(:)+1)<=rr(2) & outE(gy(:))>=rr(3)) & (outE(gy(:)+1)<=rr(4));
    C=[C,inside(:)]; %#ok<AGROW>
end,end,end,end
rect=R; overlap=O; area=A; cover=C;
end

function a=receivingAreas(e)
d=diff(e); a=kron(d',d); a=a(:);
end
function tinyHandCheck(rect,overlap,area,cover,e)
q=find(all(rect==[0 1 0 1],2),1); assert(~isempty(q));
assert(area(q)==1 && overlap(q,1)==1 && sum(overlap(q,:))==1);
assert(sum(cover(:,q))==1);
assert(receivingAreas(e)'*cover(:,q)==1);
end
function s=fieldOr(x,n,default)
if isfield(x,n),s=x.(n);else,s=default;end
end
function candidate=recheck(y,z,amp,phi,overlap,area,cover,outE,delta,rect)
C=16;L=numel(amp); level=reshape(y,L,C)'; [~,idx]=max(level,[],2); values=amp(idx)';
mass=overlap*values; selected=find(z>0); strict=1000*mass(selected)-1001*area(selected);
covered=any(cover(:,selected),2); ha=receivingAreas(outE)'*covered;
cost=sum(phi(idx));
candidate=struct('amplitude_grid',reshape(values,4,4)','amplitude_level_index',reshape(idx,4,4)', ...
 'selected_rectangles',rect(selected,:),'selected_rectangle_mass',mass(selected)', ...
 'selected_rectangle_area',area(selected)','strict_margin_numerator_1000mass_minus_1001area',strict', ...
 'halo_area',ha,'cost',cost,'ratio',ha/cost,'covered_receiving_cells',find(covered)', ...
 'recheck_ok',all(strict>0) && all(abs(strict-round(strict))<1e-9));
end
function base=axisBlockBaseline(amp,phi,overlap,area,cover,outE,delta)
best=[];
for x1=1:4,for x2=x1:4,for y1=1:4,for y2=y1:4,for l=2:numel(amp)
    v=zeros(16,1); for i=x1:x2,for j=y1:y2,v((i-1)*4+j)=amp(l);end,end
    good=(overlap*v)>(1+delta)*area; covered=any(cover(:,good),2); ha=receivingAreas(outE)'*covered; co=(x2-x1+1)*(y2-y1+1)*phi(l);
    c=struct('xblock',[x1 x2 y1 y2],'amplitude',amp(l),'halo_area',ha,'cost',co,'ratio',ha/co,'rectangle_count',sum(good));
    if isempty(best)||c.ratio>best.ratio,best=c;end
end,end,end,end,end
base=best;
end
function writeText(p,s),fid=fopen(p,'w');assert(fid>=0);fprintf(fid,'%s\n',s);fclose(fid);end
function s=notes(result)
s=sprintf(['# Correlated finite-grid probe\n\nThe MILP uses actual grid-edge rectangles, not coordinate-maximal compositions. ', ...
'Input amplitudes are rational grid values and every selected-rectangle mass, strict margin, and receiving union is independently recomputed. ', ...
'Only `exitflag=1` is an intlinprog optimum; `exitflag=2` time-limit incumbents remain candidates. ', ...
'This is a finite 4x4 family comparison and has no implication for Q, C_n, or a dimension asymptotic.\n\n', ...
'Retained rectangles: %d. Best correlated ratio: %.12g. Best axis-block constant baseline: %.12g.\n'], ...
result.rectangle_count_after_prune,result.best_candidate.ratio,result.axisblock_constant_baseline.ratio);
end
