#!/usr/bin/env python3
"""Reproduce round45 in a temporary directory; emit only a minimal receipt."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,tempfile
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]

def main():
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
    with tempfile.TemporaryDirectory(prefix='euclidean45_') as tmp:
        output=Path(tmp)/'rank.json'
        subprocess.run([sys.executable,str(HERE/'rank_probe.py'),'--output',str(output)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,env=env)
        rank=json.loads(output.read_text())
        opt=json.loads(subprocess.run([sys.executable,str(HERE/'rearrangement_probe.py')],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,env=env).stdout)
        assert rank['status']==opt['status']=='passed'
        assert rank['case_count']==19 and [b['count'] for b in rank['batches']]==[3,7,9]
        assert rank['audit']['edge_count']==806436 and rank['local_integration_tests']==225
        assert opt['case_count']==30 and opt['batch_counts']==[8,7,15]
        assert opt['max_gain_over_h_squared']=='1/6'
        assert opt['ball_benchmarks']['case_count']==11
        baseline=next(c for c in rank['cases'] if c['label']=='separated_equal_4')
        d0=next(r for r in baseline['shifted'] if r['C']=='0')
        assert F(d0['D'])==F(d0['N'])==F(1,32)
        assert F(baseline['prefix_deficit_T'])==0 and F(baseline['signed_Qhalf_minus_Cobs'])==0
        cap=[]
        for row in opt['ball_benchmarks']['cases']:
            eps=F(row['epsilon']);scaled=eps*2**32
            upper=F(-(-scaled.numerator//scaled.denominator),2**32)
            assert eps<=upper<=eps+F(1,2**32)
            cap.append(dict(n=row['n'],epsilon_upper=str(upper),epsilon_display=float(eps),
                D_over_X_certified_lower={str(C):str(max(F(0),F(1,2)-C-upper)**2) for C in (F(0),F(1,4),F(49,100))}))
        dependencies=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py']]+sorted(HERE.glob('*.py'))
        receipt=dict(status='passed',round=45,date='2026-10-07',python=sys.version.split()[0],arithmetic='Fraction exact',
            actual_rank=dict(case_count=19,dimension=1,batches=rank['batches'],audit=rank['audit'],local_integration_tests=225,
                positive_D0_count=len(rank['positive_D0_cases']),baseline={k:baseline[k] for k in ['label','X','M','Q','Cobs','prefix_deficit_T','signed_Qhalf_minus_Cobs']},baseline_D0=d0,
                max_depth=max(c['D_depth'] for c in rank['cases']),max_atoms=max(c['N'] for c in rank['cases'])),
            rearrangement=dict(case_count=30,batch_counts=opt['batch_counts'],seed=opt['seed'],max_gain_over_h_squared='1/6',scope=opt['scope']),
            high_dimension_benchmark=dict(case_count=11,batch_counts=[3,3,5],scope=opt['ball_benchmarks']['scope'],cap_certificates=cap),
            uniform_main_theorem_proved=False,sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies})
        print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
