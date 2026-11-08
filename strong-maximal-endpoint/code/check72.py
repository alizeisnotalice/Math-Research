import json,sys,re,hashlib
from pathlib import Path
from decimal import Decimal
p=Path(sys.argv[1]);py=json.load(open(p/'python/variance72.json'));wl=json.load(open(p/'wolfram/variance72.json'))
def num(s):return Decimal(re.sub(r'`[^*]*','',s).replace('*^','E'))
b={};ref=[]
for r in py['rows']:
 key=(r['N'],r['Ac'])
 if r['refinement']==1:b[key]=r
for r in py['rows']:
 if r['refinement']!=1:ref.append(abs(r['R']-b[r['N'],r['Ac']]['R']))
low={};pairs=[];cross=[]
for r in wl['rows']:
 key=(r['N'],r['Ac'])
 if r['precision']==40:low[key]=r
for r in wl['rows']:
 assert r['numeric']
 if r['precision']==80:
  key=(r['N'],r['Ac']);pairs.append(abs(num(r['R'])-num(low[key]['R'])));cross.append(abs(float(num(r['R']))-b[key]['R']))
assert wl['symbolic_hessian_identity_and_sign'] is True
assert hashlib.sha256((p/'inputs72.json').read_bytes()).hexdigest()==wl['inputs_SHA256']
assert max(ref)<2e-10 and max(cross)<2e-10 and max(pairs)<Decimal('1e-12')
assert all(r['R']>=r['positive_lower_R']-1e-10 and r['R']<r['N']+1 for r in py['rows'])
r={'pass':True,'Python_integrals':len(py['rows']),'refinement_pairs':len(ref),'refinement_max_abs':max(ref),'Wolfram_integrals':len(wl['rows']),'precision_pairs':len(pairs),'precision_max_abs':str(max(pairs)),'cross_software_pairs':len(cross),'cross_software_max_abs':max(cross),'minimum_R':min(r['R'] for r in py['rows']),'maximum_R':max(r['R'] for r in py['rows']),'scope':'numerical verification of conditional-variance identity; not an endpoint bound or interval certificate'}
Path(sys.argv[2]).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
wp=json.load(open(p/'python/window72.json'));ww=json.load(open(p/'wolfram/window72.json'));wb={(x['N'],i%4):x for i,x in enumerate(wp['rows'])};low={};err=[];diff=[]
for x in ww['rows']:
 assert x['numeric']
 if x['precision']==50:low[x['N'],x['shift_index']]=x
for x in ww['rows']:
 if x['precision']==100:
  key=(x['N'],x['shift_index']);v=num(x['tail']);err.append(abs(v-num(low[key]['tail']))/v);diff.append(abs(float(v)-wb[key]['tail'])/float(v))
assert max(err)<Decimal('1e-25') and max(diff)<2e-10
assert all(x['tail']+x['tail_prior_omission_bound']<=x['Chernoff_bound'] for x in wp['rows'])
r.update(Python_window_tail_checks=len(wp['rows']),Wolfram_window_tail_checks=len(ww['rows']),window_precision_pairs=len(err),window_precision_max_relative=str(max(err)),window_cross_pairs=len(diff),window_cross_max_relative=max(diff),maximum_prior_mass_residual=max(abs(x['mass_residual']) for x in wp['rows']))
Path(sys.argv[2]).write_text(json.dumps(r,indent=2)+'\n');print('WINDOW PASS',json.dumps({k:v for k,v in r.items() if 'window' in k or 'mass_residual' in k}))
