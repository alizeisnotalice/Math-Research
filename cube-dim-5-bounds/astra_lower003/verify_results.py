#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
import threshold_kernel as T
K=T.K;ZERO=K.ZERO;ONE=K.ONE;HERE=Path(__file__).resolve().parent
cell,z=T.cell();fine=[tuple(K.decode(t) for t in b) for b in cell['fine_cell']]
for b in K.BASE:
 assert any(fine[j][1]<=b[2*j] or b[2*j+1]<=fine[j][0] for j in range(5))
assert all(lo<z[j] and z[j]<hi for j,(lo,hi) in enumerate(fine))
assert not any(any(row['box'][2*j]==row['box'][2*j+1] for j in range(5)) for row in K.RAW)
geom=json.loads((HERE/'threshold_geometry.json').read_text());res=json.loads((HERE/'free_threshold_results.json').read_text());C=K.decode(json.loads((HERE.parent/'lower_round003/a2_equals_a3_results.json').read_text())['coefficients']);S=K.decode(geom['initial_h_slope']);assert S==C-K.decode(geom['outer_K_full'])
events=res['events'];W0=ZERO;drop=ZERO;previous=None
for e in events:
 b,a,w,t=map(K.decode,(e['beta'],e['alpha'],e['weight'],e['weighted_intercept']))
 assert ZERO<a and ZERO<w and ZERO<=b and w*b==t*a
 if previous:assert previous[0]*a<b*previous[1]
 previous=(b,a);drop+=w
 if b==ZERO:W0+=w
assert drop==K.decode(geom['total_contact_slope_drop'])
D0=K.decode(res['D_at_h0_coefficients']);assert D0==K.decode(res['B0_coefficients'])+ONE-C
assert ZERO<S and S<W0 and K.bounds(D0.a)[1]<0
assert K.decode(res['left_h_slope'])==S and K.decode(res['right_h_slope'])==S-W0
for cand in res['candidates']:
 hn,hd=map(K.decode,(cand['h_numerator'],cand['h_denominator']));assert ZERO<hd
 numerator=D0*hd+S*hn
 for e in events:
  loss=K.decode(e['weight'])*hn-K.decode(e['weighted_intercept'])*hd
  if ZERO<loss:numerator-=loss
 assert numerator==K.decode(cand['D_numerator']) and hd==K.decode(cand['D_denominator'])
 assert K.bounds(numerator.a)[1]<0
check=json.loads((HERE/'contact_check_results.json').read_text());assert check['total_area_checks']==check['total_gap_checks']==333 and check['normalized_D_numerator_matches']
checkpoint=json.loads((HERE/'direct_jet_checkpoint.json').read_text());assert checkpoint['exit_code']==1 and checkpoint['independent_direct_5D_verification_passed'] is False
out={'status':'EXACT_HOLE_AND_ALL_FIXED_REAL_H_FIRST_ORDER_CERTIFICATE_VERIFIED','fine_cell_disjointness_checked':True,'positive_initial_slope':True,'zero_event_weight_exceeds_initial_slope':True,'global_derivative_maximum_strictly_negative':True,'event_candidates_checked':len(res['candidates']),'contact_checks':333,'full5D_jet_check_complete':False,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
