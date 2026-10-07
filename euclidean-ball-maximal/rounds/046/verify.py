#!/usr/bin/env python3
"""Fresh exact round46 verification; only emit the compact certificate."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,isqrt
import hashlib,json,os,subprocess,sys,tempfile
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]

def geometry_checks():
    def exp_bounds(x,m):
        partial=sum((x**j/F(factorial(j)) for j in range(m+1)),F(0))
        next_term=x**(m+1)/F(factorial(m+1))
        return partial,partial+next_term/(1-x/F(m+2))
    assert exp_bounds(F(7,10),6)[0]>2
    assert exp_bounds(F(11,10),12)[1]<4
    assert exp_bounds(F(4),24)[1]<55
    assert 2*(110*F(11,10)+F(11,10))/F(99,100)**2<250
    batches=((10000,14400,40000),(100000,1000000,4000000),(10**8,10**10,10**12))
    rows=[]
    for batch,ns in enumerate(batches):
        for n in ns:
            source=F(3,8)-F(5,isqrt(n));observer=F(750,n+2)
            assert source>=F(13,40) and observer<F(3,40) and source-observer>F(1,4)
            assert 1/(1-F(14,5*n))<=1+F(3,n)
            assert 1/(1+F(11,10*n))<=1-F(1,n)
            assert 1-F(81,100*n)>F(99,100)
            rows.append(dict(batch=batch,n=n,source_lower=str(source),observer_upper=str(observer),gap_lower=str(source-observer),positive_transport_mass='(100*n)^(-n)/4'))
    return dict(status='passed',scope='rational parameter checks of the proved high-dimensional construction; no rare-volume sampling',case_count=9,batch_counts=[3,3,3],exponential_constant_certificates=True,cases=rows)

def main():
    with tempfile.TemporaryDirectory(prefix='euclidean46_') as folder:
        def run(name):
            out=Path(folder)/(name+'.json')
            subprocess.run([sys.executable,str(HERE/(name+'.py')),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
            result=json.loads(out.read_text());assert result['status']=='passed';return result
        search=run('prefix_probe');threshold=run('threshold_probe');geometry=geometry_checks()
        assert search['case_count']==60 and search['positive_case_count']==4
        assert [b['count'] for b in search['batches']]==[10,20,30]
        assert all(b['eligible_empty_count']==0 for b in search['batches'])
        minimal=search['minimal_witness']
        assert F(minimal['T'])==F(45,1048576) and F(minimal['T_over_X'])==F(9,14336)
        assert minimal['independent_audit']['status']==minimal['interval_certificate']['status']=='passed'
        assert threshold['batch_counts']==[3,7,9] and len(threshold['cases'])==19
        assert sum(c['threshold_interval_count'] for c in threshold['cases'])==665
        assert sum(c['pair_entries_checked'] for c in threshold['cases'])==19398
        direct=[c for c in threshold['cases'] if c['direct_pair_F'] is not None]
        assert len(direct)==3 and all(c['direct_pair_F']==c['mean_F'] for c in direct)
        dependencies=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py']]+sorted(HERE.glob('*.py'))
        receipt=dict(status='passed',round=46,date='2026-10-07',python=sys.version.split()[0],arithmetic='Fraction exact',
            search=dict(case_count=60,seed=search['seed'],batches=search['batches'],exact_checks=search['exact_checks'],
                positive_cases=[dict(label=c['label'],T_over_X=c['T_over_X'],independent_audit=c['independent_audit']) for c in search['cases'] if F(c['T'])>0]),
            minimal_witness={k:minimal[k] for k in ['atoms','alpha','X','T','T_over_X','witness','independent_audit','interval_certificate']},
            common_threshold=dict(case_count=19,batch_counts=[3,7,9],threshold_interval_count=665,pair_entries_checked=19398,
                cases=[{k:c[k] for k in ['batch','label','N','D','alpha','X','threshold_interval_count','mean_M','mean_Q','mean_Q_over_X','mean_F','mean_B','mean_diagonal','direct_pair_F']} for c in threshold['cases']],gap_tests=threshold['gap_tests']),
            high_dimension=geometry,main_theorem_proved=False,averaged_flow_budget_proved=False,
            sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies})
        print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
