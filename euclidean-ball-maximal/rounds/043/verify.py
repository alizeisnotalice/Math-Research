#!/usr/bin/env python3
"""Reproduce round43 exact ledgers without committing regenerable raw data."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,tempfile
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]

def main():
    with tempfile.TemporaryDirectory(prefix='euclidean43_') as tmpname:
        tmp=Path(tmpname)
        def run(script,name,*args):
            out=tmp/(name+'.json')
            subprocess.run([sys.executable,str(HERE/(script+'.py')),'--output',str(out),*args],check=True,
                stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
            data=json.loads(out.read_text());assert data['status']=='passed';return data
        global_data=run('global_probe','global')
        enhanced=run('global_probe','enhanced','--enhanced-only')
        geometry=run('geometry_probe','geometry')
        assert global_data['case_count']==20 and [b['count'] for b in global_data['batches']]==[4,7,9]
        assert all(v['status']=='passed' for v in global_data['independent_band_clock_verifications'])
        micro=next(c for c in global_data['cases'] if c['label']=='old41_micro_k12')
        assert F(micro['F'])==F(81,268435456) and F(micro['X'])==F(15,64)
        assert micro['coarse_cancellation_fraction']=='12/13'
        assert [c['min_S_on_last_window'] for c in enhanced['cases']]==['189/256','693/256','1449/256']
        for c in global_data['cases']:
            assert F(c['observer_square'])<=F(c['lebesgue_square'])
        scalar_keys=['X','B','F','Q','B_over_X','F_over_X','Q_over_X','coarse_positive','coarse_negative','coarse_cancellation_fraction']
        enhanced_keys=['L','observer_S_lower_bound','min_S_on_last_window','heavy_tail_volume_over_E','heavy_tail_B_over_X']
        dependencies=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/'rounds/041/square_probe.py',PROJECT/'rounds/042/row_probe.py']+sorted(HERE.glob('*.py'))
        receipt=dict(status='passed',round=43,date='2026-10-07',python=sys.version.split()[0],arithmetic='Fraction exact; rational root brackets',
            cases=20,batches=global_data['batches'],independent_verification=global_data['independent_band_clock_verifications'],
            positive_global_flow={k:micro[k] for k in scalar_keys},
            negative_rows=[dict(label=c['label'],**r) for c in global_data['cases'] for r in c['rows'] if F(r['F'])<0],
            enhanced_windows=[{k:c[k] for k in enhanced_keys} for c in enhanced['cases']],
            geometry=geometry,
            observer_square_maxima={key:dict(label=(c:=max(global_data['cases'],key=lambda c:F(c[key])))['label'],exact=c[key]) for key in ['observer_square_over_X','lebesgue_square_over_X']},
            main_bound_proved=False,uniform_global_flux_bound_proved=False,
            sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies})
        print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
