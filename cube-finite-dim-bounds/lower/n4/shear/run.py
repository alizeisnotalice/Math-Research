"""One preflight parameter then three independent fixed shear evaluations."""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
import hashlib,json,time
import engine
from validate import preflight,row_check
P=Path(__file__).resolve().parent;OLD=1.9017426111009554

def worker(t):
    row,keys=engine.evaluate(F(t));row['generator_pool_sha256']=hashlib.sha256(json.dumps(sorted(engine.old_generators()|keys),separators=(',',':')).encode()).hexdigest()
    if row['full_finite_ratio']>OLD+1e-11:(P/('candidate_'+t.replace('/','_')+'.json')).write_text(json.dumps(dict(row=row,generators=sorted(engine.old_generators()|keys)),indent=2)+'\n')
    return row

def main():
    start=time.monotonic();first=F(-1,16);validation=preflight(first);print('PREFLIGHT PASS',flush=True)
    row,keys=engine.evaluate(first);row['generator_pool_sha256']=hashlib.sha256(json.dumps(sorted(engine.old_generators()|keys),separators=(',',':')).encode()).hexdigest();validation['real_rows']=row_check(first,engine.old_generators()|keys)
    (P/'validation_results.json').write_text(json.dumps(validation,indent=2)+'\n');print('FIRST',row,flush=True);rows=[row]
    if row['full_finite_ratio']>OLD+1e-11:
        (P/'candidate_first.json').write_text(json.dumps(dict(row=row,generators=sorted(engine.old_generators()|keys)),indent=2)+'\n')
    else:
        with ProcessPoolExecutor(max_workers=3) as pool:
            for r in pool.map(worker,['-1/32','1/32','1/16']):rows.append(r);print('POINT',r,flush=True)
    out=dict(status='finite fully lattice-closed witness values only; no input/family upper bound',normalized_shear=True,rows=rows,firstpoint_preflight_passed=True,candidate=any(r['full_finite_ratio']>OLD+1e-11 for r in rows),elapsed=time.monotonic()-start,stop='candidate exact certification priority' if any(r['full_finite_ratio']>OLD+1e-11 for r in rows) else 'four fixed shear parameters complete; no optimization')
    (P/'results.json').write_text(json.dumps(out,indent=2)+'\n');print('END',out['stop'],'elapsed',out['elapsed'],flush=True)
if __name__=='__main__':main()
