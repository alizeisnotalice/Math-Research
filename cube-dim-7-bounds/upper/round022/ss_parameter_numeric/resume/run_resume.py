#!/usr/bin/env python3
"""R022 continuation: reuse a50; exactly 29 new parameters, no retries."""
from flint import arb,fmpq,ctx
from fractions import Fraction as Q
from pathlib import Path
import os,hashlib,json,time,resource,sys,datetime,traceback
D=Path(__file__).resolve().parent;P=D.parent;ctx.prec=512

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dyadic(x):
 m,e=x.man_exp();return Q(int(m)*2**int(e)) if e>=0 else Q(int(m),2**(-int(e)))
def ball(x):return arb(fmpq(x.numerator,x.denominator))
def verify():
 assert sha(P/'manifest.json')=='f091fe3e018df150a7ff0d1fc225a495e1a88e9f163e22be347ed666151f6178'
 m=json.loads((P/'manifest.json').read_text());assert all(sha(P/n)==h for n,h in m['files'].items())
 lock=json.loads((P/'sourcehash.json').read_text());old=Path(lock['old_dir'])
 assert sha(Path(lock['proof_path']))==lock['proof_sha256'] and sha(old/'manifest.json')==lock['old_manifest_sha256'] and sha(old/'result.json')==lock['old_result_sha256'] and sha(P/'run.py')==lock['run_sha256']
 rows=[json.loads(s) for s in (P/'terminal.jsonl').read_text().splitlines()];assert len(rows)==2 and rows[0]['a']==50 and rows[0]['status']=='completed' and rows[1]['a']==50 and rows[1]['status']=='failed'
 assert Q(rows[0]['lower_exact'])<Q(rows[0]['upper_exact'])
 return rows[0],old

def main():
 t=time.perf_counter();cpu=time.process_time();start=datetime.datetime.now(datetime.timezone.utc).isoformat();first,old=verify();ref=json.loads((old/'result.json').read_text());assert ref['a']==64 and ref['dimension']==7 and ref['precision_bits']==512
 oldlo=Q(ref['arb_lower_exact']);oldhi=Q(ref['arb_upper_exact']);assert oldlo<oldhi and first['old_upper_less_than_lower']==(oldhi<Q(first['lower_exact'])) and first['upper_less_than_old_lower']==(Q(first['upper_exact'])<oldlo)
 terminal=D/'terminal.jsonl';assert not terminal.exists(),'refuse repeated execution';rows=[first]
 with terminal.open('x') as out:
  for a in [v for v in range(51,81) if v!=64]:
   try:
    A=ball(Q(a+2,a));B=ball(Q(a+1,a));value=(1+A**7)*(1+B**7*(1+7*arb(a).log()))
    lo=dyadic(value.lower());hi=dyadic(value.upper());assert lo<hi
    row={'a':a,'status':'completed','lower_exact':str(lo),'upper_exact':str(hi),'display':value.str(40),'upper_less_than_old_lower':hi<oldlo,'old_upper_less_than_lower':oldhi<lo}
    out.write(json.dumps(row)+'\n');out.flush();os.fsync(out.fileno());rows.append(row)
   except Exception as exc:
    out.write(json.dumps({'a':a,'status':'failed','error':repr(exc)})+'\n');out.flush();os.fsync(out.fileno());raise
 assert len(rows)==30 and [r['a'] for r in rows]==[v for v in range(50,81) if v!=64]
 winner=min(rows,key=lambda r:Q(r['upper_exact']));whi=Q(winner['upper_exact']);assert whi<oldlo
 assert all(whi<Q(r['lower_exact']) for r in rows if r is not winner),'finite-set winner not isolated'
 grid=Q(-(-whi.numerator*10**12//whi.denominator),10**12);assert whi<grid and Q(winner['lower_exact'])>grid-Q(1,10**12)
 verify()
 return {'status':'completed','scope':'reuse a50 and evaluate 29 other registered integers, compare frozen a64; finite set only','new_evaluations':29,'reused_parameter':50,'parameters':[r['a'] for r in rows],'precision_bits':512,'old_a':64,'old_lower_exact':str(oldlo),'old_upper_exact':str(oldhi),'winner_a':winner['a'],'winner_lower_exact':winner['lower_exact'],'winner_upper_exact':winner['upper_exact'],'winner_display':winner['display'],'winner_grid_upper':str(grid),'strict_winner_over_30_others':True,'strict_improvers_vs_old_count':sum(r['upper_less_than_old_lower'] for r in rows),'source_hashes_unchanged':True,'no_global_a_optimum_claim':True,'started_utc':start,'wall_seconds':time.perf_counter()-t,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':sys.version,'executable':sys.executable,'processes':1}
if __name__=='__main__':
 try:print(json.dumps(main(),indent=2));sys.exit(0)
 except Exception as exc:traceback.print_exc(file=sys.stderr);print(json.dumps({'status':'failed','error':repr(exc)}));sys.exit(1)
