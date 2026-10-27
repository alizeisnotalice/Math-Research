#!/usr/bin/env python3
"""Check saved bounds and source hashes without repeating 512-term calculation."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys,argparse
D=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--source-root',type=Path);a=p.parse_args();m=json.loads((D/'manifest.json').read_text());assert all(sha(D/n)==v for n,v in m['files'].items());lock=json.loads((D/'sourcehash.json').read_text());root=a.source_root or Path(lock['root']);assert all(sha(root/n)==h for n,h in lock['files'].items()) and sha(D/'recheck.py')==lock['recheck_sha256']
 r=json.loads((D/'result.json').read_text());assert r['status']=='PASS' and Q(r['t'])==Q(3391749,524288) and Q(r['x_atanh'])==Q(2867461,3916037) and r['log_terms']==512 and r['sqrt_dyadic_bits']==700
 for field in ('log_bounds_outward_80_places','B_bounds_outward_80_places','F_bounds_outward_80_places'):
  l,u=map(Q,r[field]);assert 0<l<u
 lo,hi=map(Q,r['F_bounds_outward_80_places']);U=Q(r['U']);old=json.loads((root/'upper/round027/cube_cap_scalar/result.json').read_text());p=old['winner']['F']['lower'];oldlo=Q(int(p[0]))*Q(2)**int(p[1]);assert hi<U==Q(4910992105193,200000000000)<oldlo
 assert all(r[k] for k in ['B_strictly_positive','F_upper_strict_below_U','U_strict_below_R027_full_lower','lower_hash_and_ledger_match','research_paused','source_hashes_unchanged','no_optimizer_or_old_grid_replay'])
 print(json.dumps({'status':'PASS','scope':'saved outward bounds and source hashes only; no mathematical recertification rerun'}))
if __name__=='__main__':
 try:main()
 except Exception as e:print(json.dumps({'status':'FAIL','error':repr(e)}));sys.exit(1)
