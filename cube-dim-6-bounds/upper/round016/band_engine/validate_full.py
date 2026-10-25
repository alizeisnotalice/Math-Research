"""Exactly one registered oldc execution. Never called by verify_saved.py."""
import hashlib,json,subprocess,time
from pathlib import Path
from check_anchors import P,B,bounds,endpoints,contained,Q,arb,ctx
ROOT=Path('/Users/zhengzhihao/Desktop/T')
FAN=ROOT/'Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round003/direct_kernel/fan_integer_input.txt'
HA=ROOT/'Math-Research-cube6-theory/cube-dim-6-bounds/upper/round007/g_time_moments/results.json'
HT=ROOT/'Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round011/cutoff48/results.json'
X=[1000000007,998836847,892738545,807464540,467291939,96968636];DEN=1000000007
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def preflight():
 assert sha(FAN)=='7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'
 assert sha(HA)=='d77dca5f53213ca6d033700a5a658f7973602af24f6103b97b8099c3ca6d0109'
 ha=json.loads(HA.read_text());ht=json.loads(HT.read_text())
 assert ha['center_numerators']==X and ha['denominator']==DEN and ht['T']=='6/125'
 arow=next(r for r in ha['time_nodes'] if r['tau']=='1/250');trow=next(r for r in ht['rows'] if r['label']=='A007_reference')
 assert trow['point_numerators']==X and trow['denominator']==DEN
 assert arow['complete_group_terms']==trow['complete_group_terms']==967680
 assert json.loads((P/'anchors_validation.json').read_text())['status']=='PASS'
 build=json.loads((P/'build_record.json').read_text());assert sha(P/'band_engine.cpp')==build['source_sha256'];assert sha(P/'scratch/band_engine')==build['binary_sha256']
 inputs={str(p):sha(p) for p in [FAN,HA,HT,P/'registration.json',P/'anchors_validation.json']}
 (P/'inputs.json').write_text(json.dumps(dict(sha256=inputs,point_numerators=X,denominator=DEN,a='1/250',T='6/125',preflight='PASS'),indent=2)+'\n')
 return arow,trow
def check(n,arow,trow):
 assert n['status']=='all_real_nodes_supported' and n['precision_bits']==512
 assert n['full_group_terms']==967680 and n['computed_sign_pairs']==483840
 ref=arb(arow['H_arb384'])-arb(trow['H_arb384']);refd=endpoints(ref)
 assert contained(refd,n['band']),'native band must lie inside wide old Ha-HT certificate'
 nl,nu=bounds(n['band']);rl,ru=bounds(refd);assert rl<nl<=nu<ru and nl>0
 for key in ['p_a','p_T']:assert bounds(n[key])[0]>0
 elo,ehi=bounds(n['Euler_xgrad_plus5B_plusTpT_minusapa']);assert elo<=0<=ehi
 # Interval Euler independently reconstructed from saved endpoint boxes.
 def add(v,w):return v[0]+w[0],v[1]+w[1]
 def mul(q,v):return (q*v[0],q*v[1]) if q>=0 else (q*v[1],q*v[0])
 e=mul(Q(5),bounds(n['band']))
 e=add(e,mul(Q(6,125),bounds(n['p_T'])));e=add(e,mul(-Q(1,250),bounds(n['p_a'])))
 for xi,g in zip(X,n['gradient']):e=add(e,mul(Q(xi,DEN),bounds(g)))
 assert e[0]<=0<=e[1]
 return dict(status='PASS',band_old_Ha_minus_HT_reference=refd,native_band_strictly_inside_old_reference=True,full_gradient_crosschecked=False,gradient_gap='No stored full H_T48 gradient; scalar DD6 plus chain derivation and Euler only, explicit permitted validation limit.',Euler_reconstructed_fraction_endpoints=[str(v) for v in e],endpoint_densities_positive=True,new_math_bound=False)
def main():
 if (P/'execution_record.json').exists():raise SystemExit('Refuse rerun: execution record already exists')
 arow,trow=preflight();cmd=[str(P/'scratch/band_engine'),'full',str(FAN),'1','250','6','125',str(DEN),*[str(v) for v in X]]
 record=dict(attempt=1,command=cmd,source_sha256=sha(P/'band_engine.cpp'),binary_sha256=sha(P/'scratch/band_engine'),start_unix=time.time())
 (P/'execution_record.json').write_text(json.dumps(record,indent=2)+'\n')
 start=time.monotonic();r=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
 record.update(wall_seconds=time.monotonic()-start,exit_code=r.returncode,stderr=r.stderr)
 (P/'execution_record.json').write_text(json.dumps(record,indent=2)+'\n');(P/'native_result.json').write_text(r.stdout)
 if r.returncode:raise SystemExit('Native failed; stop without expansion')
 out=check(json.loads(r.stdout),arow,trow);out['source_sha256']={str(HA):sha(HA),str(HT):sha(HT)}
 (P/'full_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status='PASS',full_executions=1,kernel_seconds=json.loads(r.stdout)['kernel_seconds'],wall_seconds=record['wall_seconds'],full_gradient_crosschecked=False)))
if __name__=='__main__':main()
