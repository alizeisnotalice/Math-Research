"""Pure saved-data check: no group replay, fan determinant replay or H/G values."""
import sys
sys.dont_write_bytecode=True
import json,hashlib,math
from pathlib import Path
from fractions import Fraction as Q
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 m=json.loads((P/'manifest.json').read_text())
 for f,h in m['whitelist_sha256'].items():assert sha(P/f)==h,('local',f)
 for f,h in m['external_readonly_sha256'].items():assert sha(Path(f))==h,('external',f)
 a=json.loads((P/'adapter_results.json').read_text());assert a['all_rays_checked']==478 and a['all_simplices_checked']==13229
 assert a['distinct_images']==3840 and a['stabilizer12_exact']==12 and a['negative_orbit_representatives']==1920 and a['weight_per_representative']==24
 assert 3840*12==46080 and 1920*24==46080 and a['actual_terms']==13229*1920 and a['logical_terms']==a['actual_terms']*24
 assert not a['S2_preflight_executed'] and a['new_true_S2_pairs']==11476
 f=json.loads((P/'formal_results.json').read_text());assert len(f['configurations'])==16 and f['positive_T_evaluations']==f['full_kernel_calls']==0
 for r in f['configurations']:
  if r['configuration']%2:assert r['residual']=='0' and not r['T_substituted']
  else:
   w=list(map(Q,r['weights_exact']));c=[(Q(x),k) for x,k in r['conditions']]
   for d in range(6):assert sum(q*(0 if d<k else math.comb(d,k)*x**(d-k)) for q,(x,k) in zip(w,c))==Q(d==5)
   assert r['negative_pair_parity'] and not r['physical_G_eligibility_claim']
 abi=json.loads((P/'abi_only.json').read_text());assert abi['mode']=='ABI_ONLY_NO_KERNEL' and abi['full_kernel_calls']==abi['positive_time_evaluations']==0
 assert (abi['ulong_bytes'],abi['fmpz_bytes'],abi['arb_struct_bytes'])==(8,8,48)
 b=json.loads((P/'build_record.json').read_text());assert b['source_sha256']==sha(P/'joint_scalar_core.cpp')
 assert not (P/'full_release.json').exists()
 print('PASS saved-only: structural receipts,orbit multiplicities,16 formal identities,ABI/hash;0 full/physical T/S2 preflight/replays/writes')
if __name__=='__main__':main()
