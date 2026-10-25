"""Exact saved checkpoint chain, finite-input margin, and finite LP dual.
No search and no rerun of the original baseline or production union geometry.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, argparse
from check_free_kernel_dual import verify as verify_dual

ROOT=Path(__file__).resolve().parents[1]

def read(name):return json.loads((ROOT/name).read_text())

def run():
    names=['audit/006-satellite_n3_results.json','audit/007_results.json','audit/007b_results.json','audit/007c_results.json']
    records=[read(n) for n in names]
    intervals=[tuple(map(F,r['ratio_interval'])) for r in records]
    for lo,hi in intervals:assert lo<hi
    gains=[]
    for i in range(1,len(intervals)):
        gain=intervals[i][0]-intervals[i-1][1];assert gain>0;gains.append(str(gain))
    inp=read('lower/n3/round005c_certificate_input.json');last=records[-1]
    assert hashlib.sha256((ROOT/'lower/n3/round005c_certificate_input.json').read_bytes()).hexdigest()==last['input_sha256']=='6d77724d8c4c7eaf416c9686e9c2cb6b041dfb2e11e0e0fe769a688fa7a64bc7'
    assert len(inp['witness_keys'])==170 and last['nonempty_inner_boxes']==158 and last['nonempty_outer_boxes']==158
    assert len(last['outer_proved_empty_key_indices'])==12 and last['b']==2
    mass=sum(map(F,inp['weights']));assert mass==F(last['mass_per_index_pair'])
    finite=read('audit/007c_finite.json');P=list(map(F,inp['periods']));tau=F(finite['tau']);N=finite['N']
    loss=340*tau*sum(P[j]*P[k] for j in range(3) for k in range(j+1,3))
    calculated=(intervals[-1][0]-loss/mass)*F(2*N+1,2*N+5)**3
    assert calculated==F(finite['finite_ratio_lower']) and calculated>intervals[-2][1]
    short=F('1.7940007113851367');assert intervals[-1][0]>short
    ledger=read('results/bounds_ledger.json');assert F(ledger['tasks'][0]['interval'][0])>=intervals[-1][0]
    assert short>F(ledger['tasks'][0]['public'])
    tasks={x['id']:x for x in ledger['tasks']}
    original4=next(x for x in read('baseline/original_certificates.json') if x['dimension']==4)
    assert tasks['L4']['interval']==original4['lower_bound_interval']
    assert F(tasks['U3']['bound'])==8 and F(tasks['U4']['bound'])==16
    dual=verify_dual()
    return {'date':'2026-10-25','status':'PASS','scope':'saved exact comparisons, finite smooth margin, exact finite-grid dual; no continuous LP transfer',
      'checkpoint_files':names,'successive_exact_gain_lowers':gains,'gain_above_round4_exact_upper':str(intervals[-1][0]-intervals[0][1]),
      'simple_lower':str(short),'latest_exact_interval':last['ratio_interval'],'input_sha256':last['input_sha256'],
      'finite_smooth_ratio_lower':str(calculated),'finite_smooth_exceeds_previous_checkpoint':True,
      'new_bound_tasks':['L3'],'other_bounds_unchanged':True,'above_checked_public':True,'above_internal_round4':True,'finite_LP_dual':dual}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--save',action='store_true');a=ap.parse_args();result=run()
    if a.save:(ROOT/'results/round005_integration.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
