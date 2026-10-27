#!/usr/bin/env python3
"""Registered 30 new scalar a values, with reusable R021 exact bounds."""
from flint import arb,fmpq,ctx
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,time,resource,sys,datetime,traceback
D=Path(__file__).resolve().parent;ctx.prec=512

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dyadic(x):
 m,e=x.man_exp();return Q(int(m)*2**int(e)) if e>=0 else Q(int(m),2**(-int(e)))
def ball(x):return arb(fmpq(x.numerator,x.denominator))
def verify(lock):
 assert sha(Path(lock['proof_path']))==lock['proof_sha256']
 old=Path(lock['old_dir']);assert sha(old/'manifest.json')==lock['old_manifest_sha256'] and sha(old/'result.json')==lock['old_result_sha256']
def main():
 t=time.perf_counter();cpu=time.process_time();start=datetime.datetime.now(datetime.timezone.utc).isoformat();lock=json.loads((D/'sourcehash.json').read_text());verify(lock)
 old=json.loads((Path(lock['old_dir'])/'result.json').read_text());assert old['a']==64 and old['dimension']==7 and old['precision_bits']==512
 oldlo=Q(old['arb_lower_exact']);oldhi=Q(old['arb_upper_exact']);assert oldlo<oldhi
 terminal=D/'terminal.jsonl';assert not terminal.exists(),'refuse to overwrite or repeat batch'
 rows=[]
 with terminal.open('x') as out:
  for a in [v for v in range(50,81) if v!=64]:
   try:
    A=ball(Q(a+2,a));B=ball(Q(a+1,a)); value=(1+A**7)*(1+B**7*(1+7*arb(a).log()))
    lo=dyadic(value.lower());hi=dyadic(value.upper());assert lo<hi
    row={'a':a,'status':'completed','lower_exact':str(lo),'upper_exact':str(hi),'display':value.str(40),'upper_less_than_old_lower':hi<oldlo,'old_upper_less_than_lower':oldhi<lo}
    rows.append(row);out.write(json.dumps(row)+'\n');out.flush();os.fsync(out.fileno())
   except Exception as exc:
    out.write(json.dumps({'a':a,'status':'failed','error':repr(exc)})+'\n');out.flush();os.fsync(out.fileno());raise
 assert len(rows)==30 and [r['a'] for r in rows]==[v for v in range(50,81) if v!=64]
 winner=min(rows,key=lambda r:Q(r['upper_exact']));whi=Q(winner['upper_exact']);assert whi<oldlo
 assert all(whi<Q(x['lower_exact']) for x in rows if x is not winner),'finite winner not strictly isolated'
 grid=Q(-(-whi.numerator*10**12//whi.denominator),10**12);assert whi<grid and Q(winner['lower_exact'])>grid-Q(1,10**12)
 verify(lock)
 return {'status':'completed','scope':'30 new integer parameters, strict finite-set minimum including frozen a=64','parameters':[r['a'] for r in rows],'precision_bits':512,'old_a':64,'old_lower_exact':str(oldlo),'old_upper_exact':str(oldhi),'winner_a':winner['a'],'winner_lower_exact':winner['lower_exact'],'winner_upper_exact':winner['upper_exact'],'winner_display':winner['display'],'winner_grid_upper':str(grid),'strict_winner_over_30_others':True,'strict_improvers_vs_old_count':sum(r['upper_less_than_old_lower'] for r in rows),'source_hashes_unchanged':True,'no_global_a_optimum_claim':True,'started_utc':start,'wall_seconds':time.perf_counter()-t,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':sys.version,'executable':sys.executable,'processes':1}
if __name__=='__main__':
 try:print(json.dumps(main(),indent=2));sys.exit(0)
 except Exception as exc:traceback.print_exc(file=sys.stderr);print(json.dumps({'status':'failed','error':repr(exc)}));sys.exit(1)
