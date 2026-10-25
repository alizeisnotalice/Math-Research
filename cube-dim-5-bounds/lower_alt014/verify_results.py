from pathlib import Path
from fractions import Fraction as Q
import sys,json,hashlib
p=Path(__file__).resolve().parent;sys.path.insert(0,str(p.parent/'lower'))
from certify_rectangular import field
F,u,v,bounds,_,_=field()
def decode(terms):
 a=[Q(0)]*25
 for k,v in terms.items():i,j=map(int,k.split(','));a[5*i+j]=Q(v)
 return F(a=a)
cfile=p.parent/'lower_round003/a2_equals_a3_results.json';clo,chi=map(Q,json.loads(cfile.read_text())['ratio_interval']);rows=[]
for spec in json.loads((p/'registration.json').read_text())['inputs']:
 name=spec['name'];d=json.loads((p/(name+'_derivative.json')).read_text());inn=json.loads((p/(name+'_inner_theta_results.json')).read_text());out=json.loads((p/(name+'_outer_theta_results.json')).read_text())
 assert d['source_sha256']==hashlib.sha256((p/'first_variation.py').read_bytes()).hexdigest()
 assert inn['source_sha256']==out['source_sha256']==hashlib.sha256((p/'finite_theta.py').read_bytes()).hexdigest()
 deriv=decode(d['right_derivative_coefficients']);assert bounds(deriv.a)[1]<0
 face=[decode(x['coefficients']) for x in d['coordinate_face_contributions']]
 if name=='first_tail':
  assert face[2]==face[3]==face[4] and all(bounds(z.a)[1]<0 for z in face)
 else:
  assert face[0]==F(0) and face[1]==face[2] and face[3]==face[4] and all(bounds(z.a)[1]<0 for z in face[1:])
 assert sum((decode(x['coefficients']) for x in d['coordinate_face_contributions']),F(0))==deriv
 assert d['constant_exactly_equals_full_C'] and d['union_engines_agree'] and inn['union_engines_agree'] and out['union_engines_agree']
 assert Q(inn['ratio'])<=Q(out['ratio'])<clo<chi
 width=Q(out['ratio'])-Q(inn['ratio']);assert width<Q(3,10**17)
 for result,outer in [(inn,False),(out,True)]:
  for w in result['witnesses']:
   rr,mm=Q(w['radius']),Q(w['mass']);assert (32*rr**5>=mm) if outer else (32*rr**5<mm)
 rows.append({'name':name,'derivative_interval':list(map(str,bounds(deriv.a))),'finite_theta_ratio_interval':[inn['ratio'],out['ratio']],'exact_outer_gap_to_C_lower':str(clo-Q(out['ratio'])),'exact_enclosure_width':str(width)})
receipt={'status':'VERIFIED_TWO_NEGATIVE_EXACT_DERIVATIVES_AND_TWO_FINITE_THETA_EXCLUSIONS','rows':rows,'files':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(p.iterdir()) if f.is_file() and f.name!='completion_receipt.json'}}
(p/'completion_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(receipt['status'])
