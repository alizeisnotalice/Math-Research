% Independent endpoint optimization and explicit cross constructions.
rng(20261105,'twister');cases=2000;maxUpperExcess=0;maxAttainmentError=0;counts=zeros(1,3);
for j=1:cases
 w=.03+.97*rand;d=1.3*(1-w^2)*rand;W=1/w;
 if d<=1-w
  branch=1;L=(w+sqrt(w^2+4*d))/2;added=d;removed=0;predicted=L/sqrt(1+d);
 elseif d<1-w^2
  branch=2;L=(1-d)/w;added=L^2-w*L;removed=d-added;
  predicted=sqrt((1-d)/(2*(1-d)-w^2));
 else
  branch=3;L=w;added=0;removed=1-w^2;predicted=1;
 end
 counts(branch)=counts(branch)+1;T=W-removed/w;
 assert(T>=L-1e-12 && L>=w-1e-12 && removed>=-1e-12 && added>=-1e-12);
 mass=w*T+L^2-w*L;
 scales=[w,T,L];captured=min(T,scales).*min(w,scales) ...
  +min(L,scales).^2-min(L,scales).*min(w,scales);
 attained=max(captured./scales)/sqrt(mass);
 maxAttainmentError=max(maxAttainmentError,abs(attained-predicted));
 assert(abs((removed+added)-min(d,1-w^2))<1e-11);
 if branch<3
  ell0=(w+sqrt(w^2+4*d))/2;ell1=(1-d)/w;
  ell=unique([logspace(log10(w/2),log10(2*W),1001),w,W,ell0,ell1]);
  a=min(min(ell.^2,w*ell),1);blo=max(0,d-1+a);bhi=min(d,ell.^2-a);
  assert(all(blo<=bhi+1e-11));
  low=(a+blo)./(ell.*sqrt(1-d+2*blo));high=(a+bhi)./(ell.*sqrt(1-d+2*bhi));
  maxUpperExcess=max(maxUpperExcess,max([low,high])-predicted);
 end
end
assert(maxUpperExcess<1e-11 && maxAttainmentError<1e-11);
phaseResidual=0;phaseRows=0;maxPhaseExpansionError=0;
for A=[0,.25,1,4]
 for power=[5,7,9,11]
  t=2^-power;alpha=1-t;d=A*t^2;
  ws=alpha*sqrt(1+d)-d/(alpha*sqrt(1+d));p=-log(ws)/t;
  L=(ws+sqrt(ws^2+4*d))/2;
  phaseResidual=max(phaseResidual,abs(L/sqrt(1+d)-alpha));
  expansion=1+(1+A)*t/2+(1/3+2*A)*t^2+(1/4+3*A-A^2/4)*t^3;
  maxPhaseExpansionError=max(maxPhaseExpansionError,abs(p-expansion)/t^4);
  phaseRows=phaseRows+1;
 end
end
% Arbitrary fractional exponent: compare full b grid to the two endpoints.
generalCases=5000;maxGeneralEndpointExcess=0;
for j=1:generalCases
 gamma=.01+.98*rand;s=exp(4*rand-2);a=min(1,s)*rand;
 delta=min(.99,(1+s-2*a)*rand);blo=max(0,delta-1+a);bhi=min(delta,s-a);
 b=linspace(blo,bhi,129);f=(a+b)./(1-delta+2*b).^gamma;
 endpoints=max(f([1,end]));
 maxGeneralEndpointExcess=max(maxGeneralEndpointExcess,max(f)-endpoints);
 assert(max(f)<=endpoints+1e-12 && s^(gamma-1)*endpoints<=1+1e-12);
end
out=struct('matlab',version,'seed',20261105,'random_cases',cases,'branch_counts',counts, ...
 'scale_grid_budget',1001,'maximum_numerical_upper_excess',maxUpperExcess, ...
 'maximum_cross_attainment_error',maxAttainmentError,'phase_rows',phaseRows,'general_fractional_endpoint_cases',generalCases, ...
 'max_general_endpoint_excess',maxGeneralEndpointExcess, ...
 'maximum_phase_root_residual',phaseResidual,'maximum_phase_series_scaled_error',maxPhaseExpansionError, ...
 'scope','independent double endpoint search and finite-contact cross peak; analytic endpoint reduction proves all-input optimality');
base=fullfile(fileparts(mfilename('fullpath')),'..','results');
fid=fopen(fullfile(base,'critical_rectangle_matlab.json'),'w');fprintf(fid,'%s\n',jsonencode(out));fclose(fid);disp(out);
