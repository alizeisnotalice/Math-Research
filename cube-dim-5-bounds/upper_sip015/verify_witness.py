#!/usr/bin/env python3
"""Independent exact summation checker for saved abstract witness; no LP or solve."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,os,time
HERE=Path(__file__).resolve().parent
MOM=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/astra_upper006/verification.json')
CAP=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper-audit-cont/cube-dim-5-bounds/upper_I027/manifest.json')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();result=json.loads((HERE/'exact_witness_results.json').read_text());reg=json.loads((HERE/'registration.json').read_text())
 assert result['source_sha256']==sha(HERE/'exact_witness.py')
 assert result['moment_source_sha256']==sha(MOM)
 assert result['registration_sha256']==sha(HERE/'registration.json')
 assert result['process_results_sha256']==sha(HERE/'process_results.json')
 cap=json.loads(CAP.read_text())['caps']['1/5'];assert Q(cap)==Q(reg['Q'])
 w=result['witness'];atoms=w['atoms'];assert len(atoms)==6
 zs=[Q(x['z']) for x in atoms];weights=[Q(x['weight']) for x in atoms]
 assert len(set(zs))==6 and all(0<z<=1 for z in zs) and all(q>0 for q in weights)
 for z,row in zip(zs,atoms):assert z**20==Q(row['d'])
 m=json.loads(MOM.read_text())['moments'];targets=[Q(1)]+[Q(m[str(k)][0]) for k in range(1,5)]
 got=[sum(q*z**(20*k) for q,z in zip(weights,zs)) for k in range(5)]
 assert got==targets==list(map(Q,w['moments']))
 neg=sum(q/z**4 for q,z in zip(weights,zs));assert neg==Q(cap)==Q(w['negative_moment_Q'])
 value=sum(q/z**3 for q,z in zip(weights,zs));assert value==Q(w['objective_exact'])
 assert value>Q(3,2)>Q(17705,12500)
 # The symbolic weak-duality residual after integrating any candidate is exactly
 # alpha*(Q-neg)+sum_{k=0}^4 a_k*(m_k-got_k), hence identically zero.
 budget_identity_coefficients=[Q(cap)-neg]+[x-y for x,y in zip(targets,got)]
 assert budget_identity_coefficients==[Q(0)]*6
 source=json.loads((HERE/f"master_{w['phase']}.json").read_text())
 assert zs==[Q(x['source_node']['rational']) for x in source['dual_support']]
 out={'status':'PASS_independent_exact_witness_checker','pid':os.getpid(),'moment_residuals':[str(x-y) for x,y in zip(got,targets)],'negative_moment_residual':str(neg-Q(cap)),
 'all_six_weights_positive':True,'all_six_supports_in_zero_one':True,'objective_above_3_over_2':True,'objective_above_17705_over_12500':True,
 'budget_identity_coefficients':[str(x) for x in budget_identity_coefficients],'objective_exact':str(value),'objective_display_only':float(value),
 'witness_result_sha256':sha(HERE/'exact_witness_results.json'),'checker_source_sha256':sha(Path(__file__)),'Q_manifest_sha256':sha(CAP),
 'wall_seconds':time.monotonic()-started,'scope':'Only abstract moment-family exclusion; no real shape distribution or new I/C5 bound'}
 (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k] for k in ('status','pid','objective_display_only','wall_seconds')}))
if __name__=='__main__':main()
