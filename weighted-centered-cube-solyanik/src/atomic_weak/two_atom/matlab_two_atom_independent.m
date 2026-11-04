% Independent strict-block enumeration and weighted integration at tau=1/10.
% No clearvars or path reset.
r26_root=fileparts(mfilename('fullpath'));tau=.1;points=[-1 1];masses=[.25 15/8];
r26_blocks=[];
for i=1:2
    total=0;
    for j=i:2
        total=total+masses(j);rad=total/2;
        if rad>(points(j)-points(i))/2
            r26_blocks=[r26_blocks;points(j)-rad,points(i)+rad]; %#ok<AGROW>
        end
    end
end
r26_expected=[-9/8 -7/8;-1/16 1/16;1/16 31/16];assert(max(abs(r26_blocks-r26_expected),[],'all')<1e-13);
assert(r26_blocks(2,2)==r26_blocks(3,1)); % touching but strict endpoint excluded
r26_halo=0;
for k=1:size(r26_blocks,1),r26_halo=r26_halo+r26_weight_interval(r26_blocks(k,1),r26_blocks(k,2),tau);end
r26_cost=sum(masses.*abs(points).^(tau-1));r26_ratio=r26_halo/r26_cost;r26_tau_ratio=tau*r26_ratio;
% Independent numerical single-atom critical solve r N'-N=0.
lo=1+1e-12;hi=4;
for k=1:220
    r=(lo+hi)/2;if r26_critical(tau,r)>0,lo=r;else,hi=r;end
end
r=(lo+hi)/2;r26_single_tau=((r-1)^tau+(r+1)^tau)/(2*r);
assert(r26_tau_ratio>.87 && r26_single_tau<.865);
r26_summary=struct('version','matlab-twoatom-independent-v1-20261105','blocks',r26_blocks,'cost',r26_cost,'halo_mass',r26_halo,'ratio',r26_ratio,'tau_ratio',r26_tau_ratio,'single_numeric_tau_upper_candidate',r26_single_tau,'matlab_version',version);
r26_fid=fopen(fullfile(r26_root,'matlab_two_atom_summary.json'),'w');assert(r26_fid>=0,'Cannot write two-atom summary.');fprintf(r26_fid,'%s\n',jsonencode(r26_summary));fclose(r26_fid);
fprintf('MATLAB two atom: tauR %.15f, single candidate %.15f\n',r26_tau_ratio,r26_single_tau);
function y=r26_weight_interval(a,b,t)
if b<=0,y=((-a)^t-(-b)^t)/t;elseif a>=0,y=(b^t-a^t)/t;else,y=((-a)^t+b^t)/t;end
end
function y=r26_critical(t,r)
N=(r-1)^t+(r+1)^t;Np=t*((r-1)^(t-1)+(r+1)^(t-1));y=r*Np-N;
end
