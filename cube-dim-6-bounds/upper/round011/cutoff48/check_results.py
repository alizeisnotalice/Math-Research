#!/usr/bin/env python3
"""Independent arithmetic reconstruction from exact old G; zero kernel passes."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from time import perf_counter
import hashlib,json,os
import flint
from flint import arb,fmpz,fmpq,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def exact_g(n,d):return arb(fmpq(fmpz(n),fmpz(d)))/arb.pi()**5
start=perf_counter();r=json.loads((HERE/'results.json').read_text());m=json.loads((HERE/'input_manifest.json').read_text())
assert sha(HERE/'input_manifest.json')==r['manifest_sha256'] and sha(HERE/'run_cutoff.py')==r['source_sha256']
for inp in m['inputs']:assert sha(inp['path'])==inp['sha256'],inp['path']
assert not m['duplicate_scan']['exact_T_tag_matches']
assert r['status']=='COMPLETE_TWO_NECESSARY_CONDITIONS_ONLY_STOPPED' and r['exit_code']==0
assert r['new_H_passes']==2 and r['new_G_passes']==r['new_gradients']==r['new_spatial_points']==0
b2=json.loads(Path(m['inputs'][0]['path']).read_text());g0=json.loads(Path(m['inputs'][2]['path']).read_text())
Gs=[exact_g(b2['G']['G_times_pi5_numerator'],b2['G']['G_times_pi5_denominator']),exact_g(g0['Q_over_pi5_numerator'],g0['Q_over_pi5_denominator'])]
points=[(b2['point_numerators'],b2['denominator']),(g0['point_numerators'],g0['denominator'])]
rows=[]
for i,(row,G,(ns,d))in enumerate(zip(r['rows'],Gs,points)):
 assert row['point_numerators']==ns and row['denominator']==d
 assert row['complete_group_terms']==967680 and row['computed_sign_paired_terms']==483840
 H=arb(row['H_arb384']);K=G-H
 assert (G-arb(row['G_reconstructed_no_kernel']['ball'])).contains(0)
 assert (K-arb(row['K_T']['ball'])).contains(0)
 assert K>arb(3)/250000>arb(3)/256000
 lo=arb(12498 if i==0 else 14886)/10**9
 hi=arb(12499 if i==0 else 14887)/10**9
 assert lo<K<hi
 rows.append({'label':row['label'],'K_T':K.str(100,more=True),'strict_rational_K_interval':[f'{12498 if i==0 else 14886}/1000000000',f'{12499 if i==0 else 14887}/1000000000'],
 'margin_above_m0':(K-arb(3)/250000).str(100,more=True),'formal_point_cost':((arb(6)/125)/(64*K)).str(100,more=True),'necessary_conditions_only':True})
out={'status':'PASS_TWO_POINT_NECESSARY_CONDITIONS_ONLY','pid':os.getpid(),'elapsed_seconds':perf_counter()-start,'kernel_replayed':False,'new_kernel_passes':0,'T':'6/125','source_H_passes':2,'source_G_passes':0,'source_new_points':0,'hash_checks':'PASS','rows':rows,'global_m_T_certified':False,'new_C6_bound':False,'global_minimum_claimed':False,'python':sys.version,'python_flint':flint.__version__,'arb_bits':384,'source_sha256':sha(__file__),'results_sha256':sha(HERE/'results.json'),'input_manifest_sha256':sha(HERE/'input_manifest.json'),'exit_code':0}
(HERE/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
