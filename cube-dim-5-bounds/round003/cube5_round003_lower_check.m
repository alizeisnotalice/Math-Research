function cube5_round003_lower_check()
tic;
u=2^(1/5); v=3^(1/5);
a=[1+u^3/2, (2+u^4*v)/4+3/128, repmat((1+u)/2,1,3)];
L=[1 1 1 1 1;
1 1 1 1 2;
1 1 1 2 1;
1 1 1 2 2;
1 1 2 1 1;
1 1 2 1 2;
1 1 2 2 1;
1 1 2 2 2;
1 2 1 1 1;
1 2 1 1 2;
1 2 1 2 1;
1 2 1 2 2;
1 2 2 1 1;
1 2 2 1 2;
1 2 2 2 1;
1 2 2 2 2;
1 2 3 3 3;
1 3 2 3 3;
1 3 3 2 3;
1 3 3 3 2;
1 3 3 3 3;
2 2 2 2 2;
2 2 2 2 3;
2 2 2 3 2;
2 2 2 3 3;
2 2 3 2 2;
2 2 3 2 3;
2 2 3 3 2;
2 2 3 3 3;
2 3 2 2 2;
2 3 2 2 3;
2 3 2 3 2;
2 3 2 3 3;
2 3 3 2 2;
2 3 3 2 3;
2 3 3 3 2;
2 3 3 3 3;
2 3 4 4 4;
2 4 3 4 4;
2 4 4 3 4;
2 4 4 4 3;
2 4 4 4 4];
B=zeros(size(L,1),10);
for i=1:size(L,1)
 radius=prod(L(i,:))^(1/5)/2;
 for k=1:5
  ell=L(i,k);
  if mod(ell,2)==1
   B(i,2*k-1)=0; B(i,2*k)=min(a(k)/2,radius-(ell-1)*a(k)/2);
  else
   B(i,2*k-1)=max(0,ell*a(k)/2-radius); B(i,2*k)=a(k)/2;
  end
 end
end
ratio=32*cube5_union(B);
fprintf('Finite floating cross-check only: boxes=%d, ratio=%.17g, elapsed=%.6g s\n',size(B,1),ratio,toc);
end
function val=cube5_union(B)
if isempty(B), val=0; return; end
if size(B,2)==2
 B=sortrows(B,1); lo=B(1,1); hi=B(1,2); val=0;
 for j=2:size(B,1)
  if B(j,1)>hi, val=val+hi-lo; lo=B(j,1); hi=B(j,2); else, hi=max(hi,B(j,2)); end
 end
 val=val+hi-lo; return;
end
cuts=unique(B(:,1:2)); val=0;
for j=1:numel(cuts)-1
 lo=cuts(j); hi=cuts(j+1);
 active=B(:,1)<=lo & B(:,2)>=hi;
 val=val+(hi-lo)*cube5_union(unique(B(active,3:end),'rows'));
end
end
