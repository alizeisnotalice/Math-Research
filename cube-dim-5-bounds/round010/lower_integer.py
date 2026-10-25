"""Independent residue-block extrema; no continuous-window arrangement or area sum."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib,json,os,time
P=Path(__file__).resolve().parents[1];start=time.monotonic();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
regpath=P/'lower_round019/registration.json';resultpath=P/'lower_round019/results.json'
reg=json.loads(regpath.read_text());saved=json.loads(resultpath.read_text())
reps=[(u,v) for u,v in product(range(25),repeat=2) if (4*u+3*v)%25==0 and (-3*u+4*v)%25==0]
assert len(reps)==25
classes={};checks=0
for m in [14,15,16]:
    bywidth={}
    for width in [m,m+1]:
        rows=[]
        for a,b in product(range(25),repeat=2):
            count=sum((u-a)%25<width and (v-b)%25<width for u,v in reps)
            rows.append((count,a,b));checks+=25
        bywidth[width]={'min':min(rows)[0],'max':max(rows)[0],
                        'min_start':list(min(rows)[1:]),'max_start':list(max(rows)[1:])}
    classes[m]={'minimum_capture':bywidth[m]['min'],'maximum_capture':bywidth[m+1]['max'],
                'integer_block_checks':bywidth,'rho_open_stratum':[str(F(m,10)),str(F(m+1,10))]}
slo,shi=map(F,saved['s_interval']);ul,uu=2*slo-1,2*shi-1;assert 0<ul<uu and ul**5<2<uu**5
plans={p['N']:p for p in reg['plans']};massclasses={}
for N,p in plans.items():
    wl,wu=map(F,p['side_enclosure']);assert 0<wl<wu and wl**5<N<wu**5
    actual_l,actual_u=wl/(1+uu),wu/(1+ul)
    rl,ru=map(F,p['rho_pair']);assert rl<actual_l<actual_u<ru
    m=(10*rl).numerator//(10*rl).denominator;assert F(m)<10*rl<10*ru<F(m+1)
    massclasses[N]=m
empty=[];full=[];partial=[]
for case in saved['cases']:
    K,n,N=case['K'],case['n'],case['N'];assert K*n==N
    c=classes[massclasses[N]];mm=[c['minimum_capture'],c['maximum_capture']]
    assert case['inner_min_max']==case['outer_min_max']==mm
    if mm[1]<n:empty.append([K,n])
    elif mm[0]>=n:full.append([K,n])
    else:partial.append([K,n])
assert empty==saved['resolved_empty_cases'] and full==saved['resolved_full_cases'] and partial==saved['remaining_partial_cases']
assert checks==93750 and os.getpriority(os.PRIO_PROCESS,0)>=10 and time.monotonic()-start<1
out={'status':'PASS_independent_continuous_radius_class_extrema','normalization':'coordinates scaled by5; integer lattice residues mod25; window side L=10*rho',
     'residues':reps,'strata':classes,'actual_mass_to_stratum':massclasses,'root_endpoint_fifth_power_and_rho_brackets':True,
     'empty_case_records_verified':empty,'full_case_records_verified':full,'partial_area_cases_not_resolved_by_counts':partial,
     'scope':'Counts only; no exact partial area, union, parameter optimum or C5 improvement. Historical same-radius maxima are not claimed first discoveries.',
     'independence':'No producer atom/window/bitmask/ysweep modules or saved histogram used; direct residue block enumeration.',
     'source_sha256':{str(p.relative_to(P)):sha(p) for p in [Path(__file__),P/'round010/lower_integer_registration.md',regpath,resultpath]},
     'runtime':{'pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'membership_checks':checks,'seconds':time.monotonic()-start}}
(P/'round010/lower_integer.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','empty_case_records_verified','full_case_records_verified','runtime']},indent=2))
