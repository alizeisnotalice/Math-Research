"""Pure read: saved endpoints and provenance only, no scalar replay/writes."""
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
 dep=json.loads((P/'dependency_release.json').read_text());assert dep['status']=='A041_FROZEN_F041_PASS'
 receipt=json.loads((P/'execution_receipt.json').read_text());assert receipt['attempt']==1 and receipt['exit_code']==0 and not receipt['stderr']
 r=json.loads((P/'results.json').read_text());assert r['status']=='COMPLETE_FULL_TIME_CAP_ONLY' and len(r['configurations'])==6
 assert r['source_sha256']==sha(P/'evaluate.py') and r['kernel_calls']==0 and r['rho_values']==r['epsilon_values']==1
 assert receipt['source_sha256']==r['source_sha256'] and receipt['preparation_wait_seconds']==dep['preparation_wait_seconds']
 cs=r['configurations'];assert cs[0]['partial_fraction_residual']=='0'
 assert cs[1]['double_root_residual']==cs[1]['triple_root_residual']=='0'
 assert cs[2]['tail_derivative_residual']==cs[2]['tail_at_infinity']==cs[2]['optimizer_residual']=='0'
 assert cs[3]['rho']=='1/64' and cs[3]['A']==[992,960,928]
 assert all(bounds(v)[0]>0 for v in cs[3]['c']) and bounds(cs[3]['epsilon_star'])[0]>0
 short=bounds(cs[4]['short_norm_upper']);tail=bounds(cs[4]['tail_norm_upper']);bad=bounds(cs[5]['whole_B_norm_upper']);loose=bounds(cs[4]['loose_tail_norm_upper'])
 assert min(short[0],tail[0],bad[0])>0 and tail[1]<loose[0]
 ledger=json.loads((P/'budget_ledger.json').read_text())
 total=sum(v[1] for v in [short,tail,bad]);assert total==Q(ledger['sum_component_upper_exact'])<=Q(ledger['simple_rational_cap_norm_upper'])
 assert Q(ledger['simple_rational_cap_norm_upper'])==Q(2118205441,10**15)
 assert ledger['events']==['P:|a_v0|>=1/2','B:|a_v0|<1/2']
 assert ledger['time_ranges']==['P:(0,epsilon*)','P:[epsilon*,infinity)','B:(0,infinity)']
 assert not r['whole_face_I_certified'] and not r['g_floor_certified'] and not r['Hnorm_goal_certified'] and not r['new_C6_bound']
 print('PASS saved-only: provenance/gate,6 receipts,full-time two-event ledger,exact positive endpoint sum;0 scalar replays/0 kernels/0 writes')
if __name__=='__main__':main()
