from pathlib import Path
from fractions import Fraction as F
import json,subprocess,time,hashlib,datetime
HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'round046_gf'

def call(d,p,b,seed,n,basis=None):
 cmd=[str(HERE/'gf_check'),str(d),str(p),str(b),'2',str(seed),str(n)]+[str(x) for x in (basis or [])]
 r=subprocess.run(cmd,check=True,text=True,capture_output=True)
 return [json.loads(s) for s in r.stdout.splitlines()]

def main():
 m=json.loads((OLD/'MANIFEST.json').read_text());p=m['modulus'];K=m['K']
 pre={'small_GF8_b2':call(8,285,2,0,1,[1,6]),'GF47_b8_cost':call(47,p,8,2026102948,1),'source_C_sha256':hashlib.sha256((HERE/'gf_check.c').read_bytes()).hexdigest()}
 if not (HERE/'COST.json').exists():
  (HERE/'COST.json').write_text(json.dumps(pre,indent=2)+'\n')
 plan={'freeze_time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'written after actual small/cost runs, before any b11 bank search','d':47,'b':11,'modulus':p,'generator':2,'K':K,'seeds':[int('202610294801',16)+4096*j for j in range(4)],'max_candidates_per_bank':16,'selection':'first individually qualified basis; no rank-based reselection','source_C_sha256':pre['source_C_sha256'],'p':'99/100','n':str(10**32),'w':64,'M':40,'scales':['10/9','7/6','11/9','1',str(F(K-2,K)),str(F(K-1,K)),str(F(K+1,K)),str(F(K+2,K))]}
 if (HERE/'FROZEN_PLAN.json').exists():
  frozen=json.loads((HERE/'FROZEN_PLAN.json').read_text())
  assert all(frozen[k]==plan[k] for k in plan if k not in ['freeze_time_UTC'])
  plan=frozen
 else:
  (HERE/'FROZEN_PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
 banks=[];start=time.perf_counter()
 for j,seed in enumerate(plan['seeds']):
  attempts=call(47,p,11,seed,16);(HERE/f'BANK{j}_ATTEMPTS.json').write_text(json.dumps(attempts,indent=2)+'\n');ok=[r for r in attempts if r['passed']]
  if not ok:raise RuntimeError(f'bank {j}: all16 candidates failed')
  r=ok[0];re=call(47,p,11,0,1,r['basis'])[0]
  assert all(re[k]==v for k,v in r.items() if k not in ['attempt','seconds'])
  banks.append({'bank':j,'seed':seed,'attempts_used':len(attempts),'basis':r['basis'],'validation':r,'explicit_recheck':re})
 out={'schema':'round048.four_GF_banks.v1','d':47,'b':11,'B':2048,'q':K+1,'K':K,'modulus':p,'generator':2,'field_certificate':m['field_certificate'],'parameters':{k:plan[k] for k in ['p','n','w','M']},'banks':banks,'elapsed_seconds':time.perf_counter()-start,'FROZEN_PLAN_sha256':hashlib.sha256((HERE/'FROZEN_PLAN.json').read_bytes()).hexdigest(),'claim_scope':'Four individually qualified real GF sources; shared query, p,delta and Fourier coefficients; cross-bank ranks are measured separately.'}
 (HERE/'MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'banks':[{k:b[k] for k in ['bank','seed','attempts_used','basis']} for b in banks],'seconds':out['elapsed_seconds']},indent=2))
if __name__=='__main__':main()
