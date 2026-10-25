#!/usr/bin/env python3
"""Selected-term consistency check, not a second full pattern enumeration."""
from pathlib import Path
from itertools import permutations
from collections import Counter
from math import factorial,prod
import json,gzip,hashlib,time
from fractions import Fraction
import sympy as S
from flint import fmpq
O=Path(__file__).resolve().parent;U=O.parents[1]
tic=time.monotonic();res=json.loads((O/'results.json').read_text())
raw=gzip.decompress((O/'selected_terms.json.gz').read_bytes())
assert hashlib.sha256(raw).hexdigest()==res['uncompressed_terms_sha256']
assert hashlib.sha256((O/'selected_terms.json.gz').read_bytes()).hexdigest()==res['terms_sha256']
terms=json.loads(raw);it=iter(map(int,(U/'round003/direct_kernel/fan_integer_input.txt').read_text().split()));assert next(it)==21
cones=[]
for _ in range(21):
 de=next(it);rays=[];qs=[]
 for j in range(6):qs.append(next(it));rays.append([next(it) for _ in range(6)])
 cones.append((de,rays,qs))
count=0;rows=[]
for row,ts in zip(res['cases'],terms['cases'],strict=True):
 assert row['id']==ts['case']
 L=S.Matrix(row['normal_rows']);r=row['r'];v=S.Matrix(row['trace_direction'])
 assert int((L*L.T).det())==row['Gram_determinant']
 cols=list(zip(*row['normal_rows']));perms=sorted(set(permutations(cols)));mult=prod(factorial(k) for k in Counter(cols).values())
 assert mult==row['permutation_multiplicity'] and len(perms)==row['unique_coordinate_permutations']
 total=fmpq(0);seen=set()
 for t in ts['terms']:
  key=(t['permutation_index'],t['sign_mask'],t['cone_index']);assert key not in seen;seen.add(key)
  pc=perms[key[0]];mask=key[1];de,rays,qs=cones[key[2]]
  LS=S.Matrix([[pc[k][j]*(1 if (mask>>k)&1 else -1) for k in range(6)] for j in range(r)])
  Z=LS*S.Matrix.hstack(*[S.Matrix(ray) for ray in rays])
  J=[j for j in range(6) if any(Z[i,j] for i in range(r))];assert J==t['nonzero_normal_columns'] and len(J)==r
  M=Z[:,J];dm=int(M.det());assert dm==t['detM']
  solution=M.inv()*v;assert all(x>0 for x in solution)
  assert [int(x*dm) for x in solution]==t['cramer_numerators']
  w=fmpq(de*mult,abs(dm)*prod(qs[j] for j in range(6) if j not in J));assert str(w)==t['weight'];total+=w;count+=1
 assert str(total)==row['R_exact'] and len(seen)==row['selected_terms']
 for n in row['generic_facets']:assert sum(a*b for a,b in zip(n,row['trace_direction']))!=0
 assert sum(row['m_histogram_logical'].values())==967680
 rows.append({'id':row['id'],'selected_terms_rechecked':len(seen),'R_exact':str(total),'PASS':True})
assert count==260 and res['pattern_items']==701568 and res['logical_items']==9*967680
out={'status':'PASS_SELECTED_TERM_ARITHMETIC_AND_TOTALS','selected_terms_checked':count,'second_full_pattern_enumeration':False,'not_claimed':'Completeness is established by the producer enumeration source, its registered counts and the trace proof; this checker does not independently enumerate omitted terms.','cases':rows,'elapsed_seconds':time.monotonic()-tic,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(O/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(out['status'],count,'terms',out['elapsed_seconds'],'seconds')
