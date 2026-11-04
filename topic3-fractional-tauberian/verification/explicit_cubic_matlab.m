% Independent candidate selection; the Python ledger certifies all inequalities.
base=fullfile(fileparts(mfilename('fullpath')),'..','results');
source=jsondecode(fileread(fullfile(base,'explicit_cubic_python.json')));
rows=zeros(numel(source.profiles),6);checked=0;
for j=1:numel(source.profiles)
 p=source.profiles(j);eta=ratval(p.eta);tau=ratval(p.t_max);k=ratval(p.kappa);
 E=32*k^(2/3)+2*sqrt(k);s0=sqrt(eta)/16;
 assert(E<s0^2/4096 && k<eta^2/16 && tau<.01);
 for r=1:numel(p.inequalities)
  assert(ratval(p.inequalities(r).positive_margin)>0);checked=checked+1;
 end
 grid=linspace(0,tau,5001);
 lower=(199-200*.9+398*.9*grid+99*.9^2*grid.^2)/100;
 upper=(199-200*2+398*2*grid+99*4*grid.^2)/100;
 assert(min(lower)>0 && max(upper)<0);
 rows(j,:)=[eta,tau,k,E,s0^2/4096,max(upper)];
end
% Select dyadic eta and even binary k exponent for two prescribed t windows.
candidates=zeros(2,4);
for j=1:2
 tau=2.^(-[8,14]);tau=tau(j);
 feasible=[];
 for etaExponent=8:2:30
  eta=2^-etaExponent;
  if .1-1.8*tau-40*eta<=.09,continue;end
  if .2-2*(3*tau+40*eta)-.2*tau*(3*tau+40*eta)-.03*tau<=.1,continue;end
  exponent=fzero(@(e)32*2^(-2*e/3)+2*2^(-e/2)-eta/2^20,[40,200]);
  safeExponent=2*ceil(exponent/2);
  feasible(end+1,:)=[etaExponent,safeExponent]; %#ok<SAGROW>
 end
 [~,best]=min(feasible(:,2));candidates(j,:)=[-log2(tau),feasible(best,:),2^-feasible(best,2)];
end
out=struct('matlab',version,'positive_scalar_margins',checked, ...
 'root_polynomial_samples',numel(source.profiles)*5001,'profiles',rows,'candidate_dyadic_choices',candidates, ...
 'scope','double regressions and candidate selection; no claim of optimal mathematical constants');
fid=fopen(fullfile(base,'explicit_cubic_matlab.json'),'w');fprintf(fid,'%s\n',jsonencode(out));fclose(fid);disp(out);
function v=ratval(s)
 p=split(s,'/');v=str2double(p{1});if numel(p)==2,v=v/str2double(p{2});end
end
