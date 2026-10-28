opts=optimoptions('linprog','Display','none'); cubeRefineRows=cell(1,3);gamma=27/29;
for ii=1:3
 M=[64,256,512];M=M(ii);L=[linspace(10/9,11/9,1001)';linspace(1.30,1.40,1001)'];
 S=sin(pi*L*(1:M)).^2;t=gamma*L.^2;
 [z,val,flag]=linprog([zeros(M,1);1],[S,-t;-S,-t;ones(1,M),0],[t-1;1-t;2],[],[],zeros(M+1,1),[],opts);
 rawMinimum=min(z(1:M)); z(1:M)=max(z(1:M),0);
 d=[linspace(10/9,11/9,20001)';linspace(1.30,1.40,20001)'];
 r=struct('M',M,'rawMinimum',rawMinimum,'flag',flag,'epsilon',val,'mass',sum(z(1:M)),'amplitude',sum(sqrt(z(1:M))),'derivative',sum((1:M)'.*sqrt(z(1:M))),'dense_error',max(abs((1+sin(pi*d*(1:M)).^2*z(1:M))./(gamma*d.^2)-1)),'coefficients',z(1:M)');
 cubeRefineRows{ii}=r;fprintf('M %d eps %.10g dense %.10g amplitude %.6g derivative %.6g\n',M,val,r.dense_error,r.amplitude,r.derivative);
end
fid=fopen('refinement.json','w');fwrite(fid,jsonencode(cubeRefineRows),'char');fclose(fid);
