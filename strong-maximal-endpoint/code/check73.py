import json,sys,re,hashlib
from pathlib import Path
from decimal import Decimal
from fractions import Fraction as F
p=Path(sys.argv[1]);out=Path(sys.argv[2])
def num(s):return Decimal(re.sub(r'`[^*]*','',s).replace('*^','E'))
v={}
for stem,keys,fields in [('localize73',['N','M'],['bp','H_log','repair_log_cost']),('exponential73',['N'],['H_log']),('posterior73',['N','D_type','regime'],['actual_removed','reference_removed'])]:
 py=json.load(open(p/'python'/f'{stem}.json'));wl=json.load(open(p/'wolfram'/f'{stem}.json'));b={tuple(x[k] for k in keys):x for x in py['rows']};low={};errs=[];cross=[]
 for x in wl['rows']:
  if 'numeric' in x:assert x['numeric']
  key=tuple(x[k] for k in keys)
  if x['precision']==60:low[key]=x
 for x in wl['rows']:
  if x['precision']!=120:continue
  key=tuple(x[k] for k in keys)
  for fld in fields:
   a=num(x[fld]);bb=num(low[key][fld]);errs.append(abs(a-bb)/max(abs(a),Decimal('1e-100')));cross.append(abs(float(a)-b[key][fld])/max(abs(float(a)),1e-100))
 assert max(errs)<Decimal('1e-25') and max(cross)<2e-10
 v[stem]=dict(Python_rows=len(py['rows']),Wolfram_rows=len(wl['rows']),precision_comparisons=len(errs),precision_max_relative=str(max(errs)),cross_comparisons=len(cross),cross_max_relative=max(cross))
 if stem=='exponential73':
  for c in py['exact_certificates']:
   x=F(c['x']);target=F(c['target']);lo=sum((x**j/F(__import__('math').factorial(j)) for j in range(17)),F(0));hi=lo+x**17/F(__import__('math').factorial(17))/(1-x/18);assert str(lo)==c['lower'] and str(hi)==c['upper'];assert lo>target if c['relation']=='greater' else hi<target
  assert all(x['pass'] for x in wl['exact_certificates'])
  v[stem]['exact_rational_certificates']=4
  for x in wl['rows']:
   if x['repair_domain']:assert num(x['repair_log_cost'])>0
 if stem=='posterior73':assert hashlib.sha256((p/'python/posterior73.json').read_bytes()).hexdigest()==wl['candidate_SHA256']
v['pass']=True;v['scope']='Analytic amplitude proof plus exact rational constant inequalities; floating bounds and conditional posterior diagnostic do not certify the general sharp upper or a physical locally efficient counterexample.';out.write_text(json.dumps(v,indent=2)+'\n');print(json.dumps(v,indent=2))
