#!/usr/bin/env python3
"""Prepared one-case 128-panel auxiliary-J certificate. NOT YET RUN.

No H/R/gradient/log/exp/Gamma, norm, prior scalar or M optimization calls.
257 exact fourth-root atoms, then128 panels; atomic prefix after every atom
and panel. Explicit --resume may reuse an unfinished matching prefix only.
"""
from datetime import datetime,timezone
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import argparse
import hashlib
import json
import os
import time

HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_bytes())
def write(p,x,atomic=False):
    p=Path(p);q=p.with_suffix(p.suffix+'.tmp') if atomic else p
    with q.open('w' if atomic else 'x') as f:
        json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    if atomic:os.replace(q,p)

def root4(x,bits):
    den=1<<bits;n=x.numerator*den**4//x.denominator;k=isqrt(isqrt(n))
    lo=F(k,den);hi=lo if lo**4==x else F(k+1,den)
    assert 0<=lo and lo**4<=x<=hi**4 and hi-lo<=F(1,den)
    return lo,hi

def geometry(N):
    rows=[]
    for i in range(N):
        a,b=F(i,N),F(i+1,N);mass=(b**3-a**3)/3
        mu4=(b**7-a**7)/(7*mass);theta=(mu4-a**4)/(b**4-a**4)
        assert mass>0 and 0<=theta<=1
        rows.append((mass,mu4,theta))
    return rows

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reviewed-plan-sha',required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--resume',action='store_true');args=ap.parse_args()
    planpath=HERE/'plan.json';assert sha(planpath)==args.reviewed_plan_sha
    plan=read(planpath);assert sha(__file__)==plan['producer_sha256'] and sha(HERE/'check.py')==plan['checker_sha256']
    assert sha(HERE/'inputs.json')==plan['inputs_sha256'];inp=read(HERE/'inputs.json')
    assert plan['panels']==128 and plan['root_atoms']==257 and plan['root_bits']==192
    assert inp['epsilon']=='1/512' and inp['g']=='21/200000'
    for src in inp['source_locks'].values():assert sha(src['path'])==src['sha256']
    identity={'plan_sha256':args.reviewed_plan_sha,'inputs_sha256':plan['inputs_sha256'],
        'producer_sha256':sha(__file__),'checker_sha256':sha(HERE/'check.py'),'N':128,'bits':192}
    began=time.monotonic();N,bits=128,192
    if args.resume:
        assert sha(args.output/'inputs.json')==plan['inputs_sha256'] and sha(args.output/'plan.json')==args.reviewed_plan_sha
        state=read(args.output/'state.json');assert state['identity']==identity and state['status']!='complete'
        write(args.output/f"resume{len(state['invocations']):04d}.json",{'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),'reuses_completed_root_atoms':len(state['roots']),'reuses_completed_panels':len(state['panels'])})
        state['invocations'].append(os.getpid())
    else:
        args.output.mkdir(exist_ok=False)
        with (args.output/'inputs.json').open('xb') as f:
            f.write((HERE/'inputs.json').read_bytes());f.flush();os.fsync(f.fileno())
        with (args.output/'plan.json').open('xb') as f:
            f.write(planpath.read_bytes());f.flush();os.fsync(f.fileno())
        write(args.output/'registration.json',{'task':'upper_increment139_one_auxiliary_J_case','pid':os.getpid(),
            'utc':datetime.now(timezone.utc).isoformat(),**identity,'case_count':1,'max_root_atoms':257,
            'max_panels':128,'workers':1,'pool':None,'wall_cap':None,'new_H_R_gradient_log_exp_Gamma_calls':0,
            'old_scalar_norm_or_toy_recomputations':0,'source_locks':inp['source_locks']})
        state={'identity':identity,'status':'running_roots','roots':[],'panels':[],
            'sum_lower':'0','sum_upper':'0','invocations':[os.getpid()]}
        write(args.output/'state.json',state,True)
    eps,C,M,IH,g=map(F,(inp['epsilon'],inp['Cbar_upper'],inp['M_upper'],inp['IH_upper'],inp['g']))
    pil,piu=map(F,inp['pi_interval']);assert eps>0 and C>0 and M>0 and IH>=0 and 0<pil<=piu
    D=C/(2*eps**2);slope=M-D;geom=geometry(N)
    atoms=[('endpoint',i,D+slope*F(i,N)**4) for i in range(N+1)]
    atoms += [('mean',i,D+slope*geom[i][1]) for i in range(N)]
    try:
        # A resumed atom is only checked by powers, never recomputed with isqrt.
        for j,row in enumerate(state['roots']):
            kind,index,x=atoms[j];assert (row['kind'],row['index'],F(row['argument']))==(kind,index,x)
            lo,hi=map(F,(row['lower'],row['upper']));assert 0<=lo and lo**4<=x<=hi**4
        for j in range(len(state['roots']),len(atoms)):
            kind,index,x=atoms[j];lo,hi=root4(x,bits)
            state['roots'].append({'kind':kind,'index':index,'argument':str(x),'lower':str(lo),'upper':str(hi)})
            write(args.output/'state.json',state,True)
        state['status']='running_panels';write(args.output/'state.json',state,True)
        low,up=F(state['sum_lower']),F(state['sum_upper'])
        for i in range(len(state['panels']),N):
            mass,mu4,theta=geom[i];left,right,mean=state['roots'][i],state['roots'][i+1],state['roots'][N+1+i]
            lo=mass*((1-theta)*F(left['lower'])+theta*F(right['lower']));hi=mass*F(mean['upper'])
            assert 0<=lo<=hi;low+=lo;up+=hi
            state['panels'].append({'index':i,'mass':str(mass),'mean_z4':str(mu4),'theta':str(theta),'lower':str(lo),'upper':str(hi)})
            state['sum_lower'],state['sum_upper']=str(low),str(up);write(args.output/'state.json',state,True)
        KL,KU=5*eps**2/(2*piu),5*eps**2/(2*pil);JL,JU=KL*low,KU*up;IU=IH+JU
        result={'task':'upper_increment139','status':'one_auxiliary_integral_candidate_complete_pending_independent_primary',
            'pid':os.getpid(),'case_count':1,'N':N,'bits':bits,'root_atoms':len(state['roots']),'panels':len(state['panels']),
            'prefactor_interval':[str(KL),str(KU)],'unscaled_integral_interval':[str(low),str(up)],
            'auxiliary_J_interval':[str(JL),str(JU)],'auxiliary_J_width':str(JU-JL),'complete_I_upper':str(IU),
            'new_cost_fourth_power_upper':str(3125*IU**4/g**5),
            'below_saved135_cap_upper':JU<F(inp['saved135_cap_upper']),
            'improves_accepted_I_upper':IU<F(inp['accepted_I_upper']),
            'exact_I_upper_gain_vs_accepted':str(F(inp['accepted_I_upper'])-IU),
            'actual_J_R_or_I_lower_claim':False,'global_Cn_claim':False,
            'new_H_R_gradient_log_exp_Gamma_calls':0,'old_scalar_norm_or_toy_recomputations':0,
            'pool':None,'seconds':time.monotonic()-began,'registration_sha256':sha(args.output/'registration.json')}
        for src in inp['source_locks'].values():assert sha(src['path'])==src['sha256']
        write(args.output/'certificate.json',{'identity':identity,'roots':state['roots'],'panels':state['panels']})
        result['certificate_sha256']=sha(args.output/'certificate.json');result['all_sources_unchanged_after']=True
        write(args.output/'results.json',result);state['status']='complete';write(args.output/'state.json',state,True)
        print(json.dumps(result,sort_keys=True))
    except BaseException as exc:
        write(args.output/f"fault{len(state['invocations']):04d}.json",{'pid':os.getpid(),'exception':repr(exc),'completed_root_atoms':len(state['roots']),
            'completed_panels':len(state['panels']),'successful_prefix_preserved':True,'action':'stop;no automatic retry or refinement'})
        raise

if __name__=='__main__':main()
