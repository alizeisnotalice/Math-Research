#!/usr/bin/env python3
"""Check archived R022 bounds and provenance without recomputing Arb."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys,argparse
D=Path(__file__).resolve().parent;P=D.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--proof',type=Path);parser.add_argument('--old-dir',type=Path);a=parser.parse_args()
 m=json.loads((D/'manifest.json').read_text());assert all(sha(D/n)==h for n,h in m['files'].items())
 original=json.loads((P/'manifest.json').read_text());assert sha(P/'manifest.json')=='f091fe3e018df150a7ff0d1fc225a495e1a88e9f163e22be347ed666151f6178';assert all(sha(P/n)==h for n,h in original['files'].items())
 pre=json.loads((D/'prerunhash.json').read_text());assert all(sha(D/n)==h for n,h in pre['files'].items()) and pre['original_manifest_sha256']==sha(P/'manifest.json')
 lock=json.loads((P/'sourcehash.json').read_text());old=a.old_dir or Path(lock['old_dir']);proof=a.proof or Path(lock['proof_path']);assert sha(proof)==lock['proof_sha256'] and sha(old/'manifest.json')==lock['old_manifest_sha256'] and sha(old/'result.json')==lock['old_result_sha256']
 prior=[json.loads(x) for x in (P/'terminal.jsonl').read_text().splitlines()];assert len(prior)==2 and prior[0]['a']==50 and prior[0]['status']=='completed' and prior[1]['status']=='failed'
 new=[json.loads(x) for x in (D/'terminal.jsonl').read_text().splitlines()];assert [r['a'] for r in new]==[v for v in range(51,81) if v!=64] and all(r['status']=='completed' for r in new)
 rows=[prior[0]]+new;assert [r['a'] for r in rows]==[v for v in range(50,81) if v!=64]
 ref=json.loads((old/'result.json').read_text());lo,hi=Q(ref['arb_lower_exact']),Q(ref['arb_upper_exact']);assert lo<hi
 for r in rows:
  l,u=Q(r['lower_exact']),Q(r['upper_exact']);assert l<u and r['upper_less_than_old_lower']==(u<lo) and r['old_upper_less_than_lower']==(hi<l)
 w=next(x for x in rows if x['a']==63);wl,wu=Q(w['lower_exact']),Q(w['upper_exact']);assert wu<lo and all(wu<Q(x['lower_exact']) for x in rows if x['a']!=63)
 result=json.loads((D/'result.json').read_text());grid=Q(result['winner_grid_upper']);step=Q(1,10**12)
 assert result['status']=='completed' and result['winner_a']==63 and Q(result['winner_upper_exact'])==wu and Q(result['winner_lower_exact'])==wl and (grid/step).denominator==1 and wl>grid-step and wu<grid
 assert [r['a'] for r in rows if r['upper_less_than_old_lower']]==[62,63]
 print(json.dumps({'status':'PASS','scope':'saved exact intervals, 30 unique parameters, source hashes, strict finite-set comparisons; no Arb recomputation'}))
if __name__=='__main__':
 try:main()
 except Exception as exc:print(json.dumps({'status':'FAIL','error':repr(exc)}));sys.exit(1)
