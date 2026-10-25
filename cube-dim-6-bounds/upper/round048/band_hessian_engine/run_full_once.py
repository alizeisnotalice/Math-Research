"""One C016 oldc execution only, after the separate precheck release."""
import sys
sys.dont_write_bytecode=True
import json,subprocess,time,hashlib
from pathlib import Path
from fractions import Fraction as Q
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bounds(d):
 q=Q(2)**int(d['common_binary_exponent']);return int(d['lower_mantissa'])*q,int(d['upper_mantissa'])*q
def validate(n,old):
 assert n['full_group_terms']==old['full_group_terms']==967680 and n['computed_sign_pairs']==483840
 assert n['zero_node_terms']==old['zero_node_terms']==0 and n['confluent_terms']==old['confluent_terms']==0
 pairs=[('band',n['band'],old['band']),('p_a',n['p_a'],old['p_a']),('p_T',n['p_T'],old['p_T'])]
 pairs += [(f'gradient_{i}',v,w) for i,(v,w) in enumerate(zip(n['gradient'],old['gradient']))]
 rows=[]
 for label,v,w in pairs:
  l,u=bounds(v);a,b=bounds(w);assert max(l,a)<=min(u,b),('C016 disjoint',label)
  rows.append(dict(field=label,exact_dyadics_identical=v==w,overlap=True))
 assert len(n['hessian_upper_lex'])==21
 for v in n['hessian_upper_lex']:
  l,u=bounds(v);assert l<=u
 e=bounds(n['Euler_xgrad_plus5B_plusTpT_minusapa']);assert e[0]<=0<=e[1]
 return dict(status='PASS_IMPLEMENTATION_CROSSCHECK',existing_fields=rows,complete_Hessian_independent_reference=False,Hessian_scalar_checks=16,matching_band='[1/250,6/125]',new_Cn_bound=False,second_order_Euler_test=False)
def main():
 if (P/'full_execution.json').exists():raise SystemExit('Refuse second full execution')
 release=json.loads((P/'full_release.json').read_text());assert release['status']=='PRECHECK_ACCEPTED_FOR_ONE_OLDC_FULL'
 assert release['precheck_sha256']==sha(P/'precheck_report.md')
 reg=json.loads((P/'registration.json').read_text());build=json.loads((P/'build_record.json').read_text())
 for f,h in reg['source_runtime_input_sha256'].items():assert sha(Path(f))==h
 assert sha(P/'band_hessian.cpp')==build['source_sha256'] and sha(P/'scratch/band_hessian')==build['binary_sha256']
 assert json.loads((P/'scalar_validation.json').read_text())['status']=='PASS'
 fan=next(f for f in reg['source_runtime_input_sha256'] if f.endswith('fan_integer_input.txt'))
 cmd=[str(P/'scratch/band_hessian'),'full',fan,'1','250','6','125',str(reg['denominator']),*[str(x) for x in reg['point_numerators']]]
 rec=dict(attempt=1,command=cmd,start_unix=time.time(),preparation_and_release_wait_seconds=time.time()-P.stat().st_birthtime,source_sha256=build['source_sha256'],binary_sha256=build['binary_sha256'])
 (P/'full_execution.json').write_text(json.dumps(rec,indent=2)+'\n');t=time.monotonic()
 run=subprocess.run(cmd,capture_output=True,text=True,timeout=180);rec.update(process_wall_seconds=time.monotonic()-t,exit_code=run.returncode,stderr=run.stderr)
 (P/'full_execution.json').write_text(json.dumps(rec,indent=2)+'\n');(P/'native_hessian.json').write_text(run.stdout)
 if run.returncode:raise SystemExit('Whole invocation failed;stop,no retry')
 n=json.loads(run.stdout);old=json.loads((P.parents[1]/'round016/band_engine/native_result.json').read_text());v=validate(n,old)
 (P/'full_validation.json').write_text(json.dumps(v,indent=2)+'\n');print(json.dumps(dict(status=v['status'],kernel_seconds=n['kernel_seconds'],wall_seconds=rec['process_wall_seconds'],existing_9_fields_exact=all(x['exact_dyadics_identical'] for x in v['existing_fields']),Hessian_independent_full_reference=False)))
if __name__=='__main__':main()
