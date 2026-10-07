#!/usr/bin/env python3
"""Fresh temporary reproduction of the round44 truncation and slice certificates."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,tempfile
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]

def main():
    with tempfile.TemporaryDirectory(prefix='euclidean44_') as folder:
        def run(name):
            out=Path(folder)/(name+'.json')
            subprocess.run([sys.executable,str(HERE/(name+'.py')),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
            data=json.loads(out.read_text());assert data['status']=='passed';return data
        tails=run('tail_probe');slices=run('slice_probe')
        assert tails['case_count']==19 and slices['case_count']==28
        assert [b['count'] for b in tails['batch_maxima']]==[3,7,9]
        endpoints=sum(c['prefix_Lipschitz_audit']['edge_endpoint_count'] for c in tails['cases'])
        assert endpoints==1612872
        best=max(tails['cases'],key=lambda c:F(c['prefix_Lipschitz_audit']['max_difference']))
        assert F(best['prefix_Lipschitz_audit']['max_difference'])==F(511848731,4076863488)
        last=next(c for c in tails['cases'] if c['label']=='old41_chain_L96')
        keys=['t','Qtail_over_X','t_Qtail_over_X','H_over_Qtail','n1_contraction_factor','rho_sup']
        dependencies=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py']]+sorted(HERE.glob('*.py'))
        # Retain only one extremal case per batch/threshold, not the regenerable full ledgers.
        maxima=[]
        for batch in tails['batch_maxima']:
            maxima.append(dict(batch=batch['batch'],count=batch['count'],thresholds=[dict(t=r['t'],active_count=r['active_count'],maximum={k:r['maximum'][k] for k in ['H_over_Qtail','t_Qtail_over_X','rho_sup']}) for r in batch['thresholds']]))
        receipt=dict(status='passed',round=44,date='2026-10-07',python=sys.version.split()[0],arithmetic='Fraction exact',
            case_count=19,batch_maxima=maxima,prefix_audit=dict(edge_endpoint_count=endpoints,max_label=best['label'],max_difference=best['prefix_Lipschitz_audit']['max_difference'],proved_n1_bound='1/2'),
            shifted_moment_and_stop_loss_checks=19,deep_case=dict(label=last['label'],N=last['N'],D=last['D'],max_source_A=last['max_source_A'],thresholds=[{k:r[k] for k in keys} for r in last['thresholds']],rank_budget=last['rank_budget']),
            slice_checks=slices,main_bound_proved=False,dimension_independent_shift_proved=False,
            sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies})
        print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
