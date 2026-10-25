"""Read only; no Arb configuration replay, kernels, enumeration, or writes."""
import sys
sys.dont_write_bytecode=True
import json,hashlib
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bounds(v):
 e=Q(2)**int(v['common_binary_exponent']);return int(v['lower_mantissa'])*e,int(v['upper_mantissa'])*e
def main():
 m=json.loads((P/'manifest.json').read_text())
 for f,h in m['whitelist_sha256'].items():assert sha(P/f)==h,('local',f)
 for f,h in m['external_readonly_sha256'].items():assert sha(Path(f))==h,('external',f)
 reg=json.loads((P/'registration.json').read_text());res=json.loads((P/'results.json').read_text());a=json.loads((P/'audit_dependency.json').read_text())
 assert a['status']=='F035_PASS_DEPENDENCY_DISCHARGED' and a['extra_configurations']==0
 assert res['source_sha256']==sha(P/'check_diagnostic.py')
 assert len(res['configurations'])==reg['planned_configurations']==7<=12
 assert all(res[k]==0 for k in ['kernel_calls','old_enumerations','grids','T_searches','MATLAB_calls'])
 assert 192+480+320==992 and 2+990==992
 for c in res['configurations'][:3]:
  A,C=Q(c['A_times_pi6']),Q(c['C_times_pi6'])
  assert C-960*A==Q(c['C_minus960A_times_pi6'])>0
  assert bounds(c['M'])[0]>0
 for c in res['configurations'][3:6]:
  u=Q(c['u']);r=Q(c['r_squared']);dv=Q(c['d_v_squared']);dw=Q(c['d_w_squared'])
  assert 1/r*(1/dv+1/dw)==Q(c['leading_coefficient_over_M3_exact'])==Q(4,9)/(1-u*u)**2
  assert bounds(c['collapsed_over_resolved'])[0]>700
 c=res['configurations'][6]
 assert Q(c['delta'])==Q(1,64) and c['ordered_pair_count']==992 and c['other_ordered_pairs']==990
 assert bounds(c['line_gap_a'])[0]>0 and bounds(c['other_plane_gap_g'])[0]>0
 assert bounds(c['total_L6over5_upper'])[1]<Q(1,11250)
 assert c['target_ordered_pairs']==2
 print('PASS saved-only: hashes,F035 dependency,7 receipts,exact identities,all992 pair accounting,single-tube endpoints;0 scalar replays/0 kernels/0 writes')
if __name__=='__main__':main()
