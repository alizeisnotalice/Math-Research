% MATLAB R2026a independent outer-weight local blowup cross-check.
% Outer halo depends ONLY on Lebesgue lengths. Weight primitives are used later
% to measure the same H/E; no w appears in the halo predicate.
root = fileparts(mfilename('fullpath'));
taus = [0.15 0.25 0.5 0.75 0.9]; records = struct([]); count=0;
for caseId=1:24
    tau=[taus(mod(caseId-1,5)+1), taus(mod(2*caseId,5)+1)];
    if caseId>18, z=[-0.10-0.01*mod(caseId,3),0.11+0.015*mod(caseId,3)];
    else, z=[-0.72+0.06*mod(caseId,4),0.40+0.08*mod(3*caseId,4)]; end
    c0=.2+.05*mod(caseId,3); c=[.35+.11*mod(caseId,5),.55+.13*mod(2*caseId,5)];
    sep=abs(z(1)-z(2)); alpha=.9;
    for j=1:2
        for k=4:9
            epsv=2^(-k); if 5*epsv>=sep, continue; end
            E=[z(j)-1.5*epsv,z(j)-.05*epsv; z(j)+.05*epsv,z(j)+1.5*epsv];
            H=outerHaloFloat(E,alpha);
            pure=intervalMass(H,0,c,z,tau,j)/intervalMass(E,0,c,z,tau,j);
            mixed=intervalMass(H,c0,c,z,tau,0:2)/intervalMass(E,c0,c,z,tau,0:2);
            componentRatios=zeros(1,3);
            for q=0:2, componentRatios(q+1)=intervalMass(H,c0,c,z,tau,q)/intervalMass(E,c0,c,z,tau,q); end
            count=count+1; records(count).caseId=caseId; records(count).component=j;
            records(count).k=k;records(count).epsilon=epsv;records(count).tau=tau(j);
            records(count).mixtureRatio=mixed;records(count).pureRatio=pure;
            records(count).difference=mixed-pure;
            records(count).scaledAbsDifference=abs(mixed-pure)/epsv^(1-tau(j));
            records(count).convexPass=mixed<=max(componentRatios)+1e-11;
            records(count).haloIntervals=size(H,1);
        end
    end
end
out.matlabVersion=version;out.examples=count;
out.definition='Outer H_alpha(E): Lebesgue centered intervals; signed-power weights only measure H/E';
out.records=records;out.allFinite=all(isfinite([records.mixtureRatio]))&&all(isfinite([records.pureRatio]));
out.allConvexPass=all([records.convexPass]);out.maxScaledAbsDifference=max([records.scaledAbsDifference]);
fid=fopen(fullfile(root,'matlab_outer_summary.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);
fprintf('MATLAB %s outer examples=%d finite=%d convex=%d maxScaled=%.6g\n',version,count,out.allFinite,out.allConvexPass,out.maxScaledAbsDifference);

function value=signedPower(x,t), value=sign(x).*abs(x).^t; end
function value=primitive(a,b,c0,c,z,tau,indices)
 value=0;if any(indices==0),value=value+c0*(b-a);end
 for q=1:numel(c),if any(indices==q),value=value+c(q)*(signedPower(b-z(q),tau(q))-signedPower(a-z(q),tau(q)))/tau(q);end,end
end
function value=intervalMass(I,c0,c,z,tau,indices)
 value=0;for q=1:size(I,1),value=value+primitive(I(q,1),I(q,2),c0,c,z,tau,indices);end
end
function out=mergeIntervals(I)
 I=sortrows(I);out=[];for q=1:size(I,1),a=I(q,1);b=I(q,2);if b<=a,continue;end
  if ~isempty(out)&&a<=out(end,2),out(end,2)=max(out(end,2),b);else,out=[out;a b];end
 end
end
function H=outerHaloFloat(E,alpha)
 ends=E';ends=ends(:)';ell=sum(E(:,2)-E(:,1));lo=ends(1)-ell/(2*alpha);hi=ends(end)+ell/(2*alpha);pieces=E;
 for e=ends
  cuts=unique(sort([lo hi (e+ends)/2]));cuts=cuts(cuts>=lo&cuts<=hi);
  for q=1:numel(cuts)-1
   a=cuts(q);b=cuts(q+1);ga=G(a,e,E,alpha);gb=G(b,e,E,alpha);tol=1e-12;
   if abs(ga)<tol,ga=0;end;if abs(gb)<tol,gb=0;end
   if ga>0&&gb>0,pieces=[pieces;a b];elseif ga>0,pieces=[pieces;a a+(b-a)*ga/(ga-gb)];elseif gb>0,pieces=[pieces;a+(b-a)*(-ga)/(gb-ga) b];end
  end
 end;H=mergeIntervals(pieces);
end
function y=G(x,e,E,alpha)
 a=min(e,2*x-e);b=max(e,2*x-e);m=0;for q=1:size(E,1),m=m+max(0,min(b,E(q,2))-max(a,E(q,1)));end;y=m-2*alpha*abs(x-e);
end
