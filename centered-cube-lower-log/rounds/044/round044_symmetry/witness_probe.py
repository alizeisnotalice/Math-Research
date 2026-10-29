"""Same frozen tests: distinguish short/long-window evidence, no new samples."""
from pathlib import Path
import json,time
import numpy as np
import probe as p
HERE=Path(__file__).resolve().parent

def main():
 target=HERE/'WITNESS.json';assert not target.exists();start=time.time();saved=json.loads((HERE/'SUMMARY.json').read_text());rows=[]
 lengths=np.array(p.C['lengths']);short=lengths<1;long=~short
 for case in saved['cases']:
  n=case['n'];frozen=json.loads((HERE/f'FROZEN_n{n}.json').read_text());assert p.sha(HERE/f'FROZEN_n{n}.json')==case['frozen_SHA256']
  s=np.array(frozen['selector']);lam=frozen['lambda_'];rng=np.random.default_rng(p.C['test_seed']+n)
  hist=np.zeros(len(lengths),dtype=int);counts=np.zeros(4,dtype=int)
  for at in range(0,p.C['test_N'],p.C['batch']):
   x=rng.random((min(p.C['batch'],p.C['test_N']-at),n));feat,_,_=p.features(x);v=p.old.value(feat,s)
   hit=v.max(-1)>lam+p.C['guard'];sh=v[:,short].max(-1)>lam+p.C['guard'];lo=v[:,long].max(-1)>lam+p.C['guard']
   hist+=np.bincount(v.argmax(-1)[hit],minlength=len(lengths));counts+=np.array([hit.sum(),sh.sum(),lo.sum(),(sh&~lo).sum()])
  assert counts[0]==case['counts'][0] and hist.sum()==counts[0]
  rows.append(dict(n=n,lambda_=lam,all_gate_success=int(counts[0]),any_short_success=int(counts[1]),any_long_success=int(counts[2]),short_only_success=int(counts[3]),successful_argmax_histogram=hist.tolist(),frozen_SHA256=case['frozen_SHA256']))
 out=dict(cases=rows,lengths=lengths.tolist(),script_SHA256=p.sha(__file__),source_SHA256=p.sha(p.__file__),seconds=time.time()-start,
 scope='Same frozen FP query sequence; scope diagnostic only, not interval certification or all-threshold/global scale claim')
 p.dump(target,out);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
