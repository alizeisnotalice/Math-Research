% Frozen LP search: nonnegative spectrum and common variance target.
M=128; gamma=27/29; results=struct([]); counter=0;
windows={ [10/9,11/9], [10/9,11/9;1.30,1.40], [10/9,11/9;1.50,1.60], [10/9,11/9;16/9,17/9] };
opts=optimoptions('linprog','Display','none');
for caseid=1:numel(windows)
 W=windows{caseid}; L=[];
 for j=1:size(W,1),L=[L;linspace(W(j,1),W(j,2),101)'];end
 S=sin(pi*L*(1:M)).^2; target=gamma*L.^2;
 for budget=[2,8,32]
  % |1+S*b-target| <= epsilon*target, sum b <= budget.
  A=[S,-target;-S,-target;ones(1,M),0]; rhs=[target-1;1-target;budget];
  [z,obj,flag]=linprog([zeros(M,1);1],A,rhs,[],[],zeros(M+1,1),[],opts);
  counter=counter+1; results(counter).caseid=caseid;results(counter).budget=budget;results(counter).flag=flag;
  if flag>0
   dense=[];for j=1:size(W,1),dense=[dense;linspace(W(j,1),W(j,2),10001)'];end
   err=abs((1+sin(pi*dense*(1:M)).^2*z(1:M))./(gamma*dense.^2)-1);
   results(counter).epsilon=obj;results(counter).dense_error=max(err);results(counter).mass=sum(z(1:M));
   results(counter).amplitude_proxy=sum(sqrt(z(1:M)));results(counter).derivative_proxy=sum((1:M)'.*sqrt(z(1:M)));
   results(counter).coefficients=z(1:M)';
  end
 end
end
fid=fopen('results.json','w');fwrite(fid,jsonencode(results),'char');fclose(fid);
for k=1:numel(results),r=results(k);fprintf('case %d budget %g flag %d eps %.9g dense %.9g mass %.9g\n',r.caseid,r.budget,r.flag,r.epsilon,r.dense_error,r.mass);end
