#!/usr/bin/env python3
"""Read-only source and saved-result consistency; no Arb recomputation."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys,argparse
D=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--proof',type=Path);a=p.parse_args()
 m=json.loads((D/'manifest.json').read_text());assert all(sha(D/n)==h for n,h in m['files'].items())
 s=json.loads((D/'sourcehash.json').read_text());assert sha(a.proof or Path(s['proof_path']))==s['proof_sha256'] and sha(D/'run.py')==s['run_sha256']
 r=json.loads((D/'result.json').read_text());assert (r['a'],r['dimension'],r['precision_bits'])==(64,7,512)
 lo,hi,grid=map(Q,(r['arb_lower_exact'],r['arb_upper_exact'],r['grid_upper']));unit=Q(1,10**12)
 assert lo<hi<grid and lo>grid-unit and (grid/unit).denominator==1
 cubic=1+Q(7,10)+Q(7,10)**2/2+Q(7,10)**3/6
 upper=(1+Q(33,32)**7)*(1+Q(65,64)**7*(1+42*Q(7,10)))
 assert cubic==Q(r['cubic_exp_lower_exact'])>2 and upper==Q(r['rational_upper_using_log_bound'])
 assert 1+Q(33,32)**7<Q(57,25) and Q(65,64)**7<Q(64,57)
 assert hi<upper<Q(10013,125)<81<128
 print(json.dumps({'status':'PASS','scope':'saved endpoints, exact rational inequalities and hashes; not an independent Arb run or proof audit'}))
if __name__=='__main__':
 try:main()
 except Exception as e:print(json.dumps({'status':'FAIL','error':repr(e)}));sys.exit(1)
