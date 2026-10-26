#!/usr/bin/env python3
"""Only after review: six new scalar certificate objects; no H kernel imports."""
import sys
sys.dont_write_bytecode=True
import argparse,hashlib,json,os,time
from pathlib import Path
from fractions import Fraction as F
import arith as A

HERE=Path(__file__).resolve().parent
SCRATCH=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_join132')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,x):
    p=Path(p);q=Path(str(p)+'.tmp')
    with q.open('w') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    q.replace(p)
def locks(reg):
    for p,h in {**reg['input_locks'],**reg['code_locks'],**reg['protected123_locks']}.items():assert sha(p)==h,('source changed',p)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reviewed-registration-sha',required=True);ap.add_argument('--primary-authorization',required=True);a=ap.parse_args()
    assert a.primary_authorization=='upperjoin132-finite-budget-approved'
    rp=HERE/'registration.json';assert sha(rp)==a.reviewed_registration_sha;reg=read(rp);locks(reg);assert sha(HERE/'inputs.json')==reg['inputs_sha256']
    inp=read(HERE/'inputs.json');eps=F(inp['epsilon']);q=F(5,4);UH=F(inp['UH']);EU=F(inp['E_U']);C3=F(inp['C3_U']);C1=F(inp['C1_U']);piL=F(inp['pi_L']);P5=F(inp['P3_U'])/(64*piL)
    assert eps>0 and C3>0 and C1>0 and P5>0
    path=SCRATCH/'results.json';began=time.monotonic();out=read(path) if path.exists() else dict(status='finite_prefix',registration_sha256=a.reviewed_registration_sha,cases={},pid=os.getpid())
    assert out['registration_sha256']==a.reviewed_registration_sha and set(out['cases'])<=set(reg['finite_cases'])
    assert out['status']!='complete','successful whole finite batch must not be repeated'
    # Shared new scalars are durable and reused on a later finite-prefix resume.
    if 'common' in out:
        saved=out['common'];Epower=(F(saved['Epower']['lower']),F(saved['Epower']['upper']))
        mass=(F(saved['R_mass_upper_formula']['lower']),F(saved['R_mass_upper_formula']['upper']))
        assert saved['P5_upper_constant']==str(P5) and saved['UH']==str(UH) and saved['g']==inp['g']
    else:
        Epower=A.root(EU**5,4);mass=A.div(A.log(1+256*eps**2),512*piL)
    base=A.add(UH,Epower)
    common=dict(Epower=A.pack(Epower),R_mass_upper_formula=A.pack(mass),P5_upper_constant=str(P5),UH=str(UH),g=inp['g'])
    def check(name,ball,extra,valid_upper=True):
        lo,hi=ball;g=F(inp['g']);Tlo,Thi=F(inp['threshold_L']),F(inp['threshold_U']);assert 0<lo<=hi
        out['cases'][name]=dict(certificate_object=A.pack(ball),is_actual_I_upper=valid_upper,exact32_pass=valid_upper and 3125*hi**4<32**4*g**5,
            certificate_strictly_above_threshold=lo>Thi,new_upper_below_accepted_exact_I=hi<F(inp['existing_exact_I_upper']),
            new_upper_below_accepted_I_cap=hi<F(inp['existing_I_cap']),S_interval=dict(lower=str(Tlo-hi),upper=str(Thi-lo)),**extra)
        out.update(common=common,pid=os.getpid(),seconds=time.monotonic()-began);write(path,out)
        print(json.dumps(dict(case=name,upper_display=float(hi),exact32_pass=out['cases'][name]['exact32_pass'])),flush=True)
    if 'zero_cross_floor' not in out['cases']:check('zero_cross_floor',base,dict(scope='lower floor of this fixed-UH/Epower upper-certificate family only;not actualI lower bound or validI upper'),valid_upper=False)
    def B(t):
        assert t[0]>=eps
        term=A.mul(C3/2,A.sub(eps**-2,A.div(1,A.power(t,2))))
        return A.add(A.add(term,A.mul(C1,A.log(A.div(t,eps)))),A.div(P5,A.mul(4,A.power(t,4))))
    todo=[('T_'+v.replace('/','over'),A.interval(v)) for v in reg['rational_Ts']]
    if 'analytic_stationary' not in out['cases']:
        discr=A.root(C3*C3+4*C1*P5,2);z=A.div(2*P5,A.add(C3,discr));stationary=A.root(z,2)
        if stationary[1]<=eps:stationary=A.interval(eps)
        elif stationary[0]<eps:raise AssertionError('stationary straddles eps;explicit finer authorization needed')
        todo.append(('analytic_stationary',stationary))
    for name,t in todo:
        if name in out['cases']:continue
        b=B(t);assert b[0]>0;combined=A.add(base,A.mul(q,A.mul(A.root(b,4),mass)))
        check(name,combined,dict(T=A.pack(t),H_sup_upper_formula=A.pack(b),mixed_upper=A.pack(A.mul(q,A.mul(A.root(b,4),mass))),scope='new nonlinear upper certificate with global positiveface sup and traceL1'))
    assert set(out['cases'])==set(reg['finite_cases']);locks(reg)
    out.update(status='complete',all_six_new_cases=True,new_H=0,new_R_pointwise=0,new_gradient=0,old_baseline_recomputed=False,global_C5_claim=False,wall_seconds=time.monotonic()-began,
        fixed_E_zero_cross_floor_excluded=out['cases']['zero_cross_floor']['certificate_strictly_above_threshold'],
        continuous_sup_certificate_family_excluded=out['cases']['analytic_stationary']['certificate_strictly_above_threshold'],
        certificate_family_scope='fixed112_UH,076_E_U,global076face_sup+057P5_tail,g105;never actualI/model/C5 exclusion')
    write(path,out);print(json.dumps(dict(status='complete',cases=6,results_sha256=sha(path),seconds=out['wall_seconds'])),flush=True)

if __name__=='__main__':main()
