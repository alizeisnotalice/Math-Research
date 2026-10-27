#!/usr/bin/env python3
"""Saved-result accounting and exact witness checks; no kernel point calls."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    reg=json.loads((ROOT/'registration.json').read_text())
    for p,h in reg['sources_sha256'].items(): assert sha(ROOT/p)==h
    d=json.loads((ROOT/'diagnostics.json').read_text()); stages=[json.loads((ROOT/f'{p}-receipt.json').read_text()) for p in ('first8','rest56')]
    starts=[v for s in stages for v in s['starts']]
    assert sorted(v['start_id'] for v in starts)==list(range(64))
    assert len({v['handle'] for v in starts})==64
    calls=0
    for v in starts:
        full=json.loads((ROOT/'starts'/f"{v['start_id']:03d}.json").read_text())
        assert full['status']=='finished' and full['full_point_calls']==len(full['trace'])<=256
        assert full['solver_nfev']==full['solver_njev']==full['full_point_calls']
        assert all('error' not in row for row in full['trace'])
        calls+=full['full_point_calls']
    assert calls==d['optimization_calls']==868
    assert d['diagnostic_full_point_calls']==8*27==216
    assert d['high_precision_same_point_subset_verification_calls']==8*3==24
    assert d['total_budget_charged_calls']==1108<=32768
    assert d['total_full_point_calls']==1084 and not d['new_shapes']
    vv=[[1]+[1 if (b>>k)&1 else -1 for k in range(6)] for b in range(64)]
    audits=[]
    def det(a): return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
    for c in d['candidate_bundles']:
        p=c['rounded_rational_point']; x=list(map(F,p['rational'])); witness=p['S3']; js=witness['body_indices']; co=list(map(F,witness['coefficients']))
        assert all(sum(co[k]*vv[js[k]][i] for k in range(3))==x[i] for i in range(7))
        orig=[F(float.fromhex(h)) for h in c['original_optimized_point']['point']['floathex']]
        # Verify exact floathex optimizer point in the same body span as its rounded point.
        import itertools
        for rows in itertools.combinations(range(7),3):
            a=[[vv[j][r] for j in js] for r in rows]; dd=det(a)
            if dd: break
        coeff=[]
        for k in range(3):
            ac=[v.copy() for v in a]
            for rr in range(3): ac[rr][k]=orig[rows[rr]]
            coeff.append(det(ac)/dd)
        assert all(sum(coeff[k]*vv[js[k]][i] for k in range(3))==orig[i] for i in range(7))
        near=list(map(F,c['fixed_nearby_offS3_probe']['rational']))
        assert len({abs(v) for v in near})==7
        assert len(p['blocks'])==8 and p['conditioning']['all_shapes_retained']==8192
        assert all(len(c[key]['same_point_high_condition_80digit'])==7 for key in ('original_optimized_point','rounded_rational_point','fixed_nearby_offS3_probe'))
        audits.append({'handle':c['source_start_handle'],'original_exact_floathex_point_in_S3':True,'original_floathex_as_exact_rationals':[str(z) for z in orig],'body_indices':js,'original_exact_coefficients':[str(z) for z in coeff],'rounded_point_witness_verified':True,'fixed_near_point_distinct_abs':7,'fixed_near_point_off_S3':True})
    out={'status':'pass','no_kernel_evaluations':True,'original_and_rounded_8_candidates_all_in_S3':True,'nearby_8_candidates_all_off_S3':True,'counted_total_budget_calls':1108,'original_cache_shape_count':8192,'all_start_handles_terminal':True,'new_C7_bound':False,'audits':audits}
    (ROOT/'verification.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps({k:v for k,v in out.items() if k!='audits'},indent=2))
if __name__=='__main__': main()
