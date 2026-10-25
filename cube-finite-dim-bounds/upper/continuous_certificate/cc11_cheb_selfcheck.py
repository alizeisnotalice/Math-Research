"""Independent finite algebra/consistency checks; not a substitute for external audit."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import runpy,json,time,hashlib
from flint import arb,arb_mat,ctx
OUT=Path(__file__).resolve().parent;C=runpy.run_path(str(OUT/'cc11_cheb_check.py'));ctx.prec=384
N=C['N'];n=N-1;start=time.perf_counter();count=0
# Exact polynomial identity: degree <=20, verified at 21 distinct rational points.
for l in range(N):
 u=F(l,n);x=2*u-1;T=[F(1),x]
 for k in range(2,N):T.append(2*x*T[-1]-T[-2])
 for k in range(N):
  val=sum((C['BF'][j][k]*comb(n,j)*u**j*(1-u)**(n-j) for j in range(N)),F(0))
  assert val==T[k];count+=1
# DCT at exact Chebyshev nodes must reproduce the basis matrix.
Q=arb_mat([[(arb.pi()*k*(2*j+1)/(2*N)).cos() for k in range(N)] for j in range(N)])
D=C['A']*Q
for i in range(N):
 for j in range(N):assert D[i,j].contains(int(i==j)) and D[i,j].rad()<arb(2)**-300
# End-to-end off-node checks at rational u=1/3,v=2/3 on the four representative patches.
rows=[]
for i,j in [(15,7),(15,15),(8,3),(0,0)]:
 V,hit=C['node_matrix'](i,j);B=C['TRANSFORM']*V*C['TRANSFORM'].transpose()
 u=F(1,3);v=F(2,3)
 p=sum((B[a,b]*C['aq'](F(comb(n,a))*u**a*(1-u)**(n-a)*F(comb(n,b))*v**b*(1-v)**(n-b)) for a in range(N) for b in range(N)),arb(0))
 exact=C['ns']['point_divided_difference'](F(i,16)+u/16,F(j,16)+v/16)
 assert abs(p-exact)<C['INTERPOLATION_ERROR']
 rows.append({'patch':[i,j],'polynomial':str(p),'direct_H':str(exact),'absolute_difference_upper':str(abs(p-exact).upper()),'less_than_uniform_error':True})
result=json.loads((OUT/'cc11_cheb_all_result.json').read_text())
assert result['coverage_complete'] and result['config']==C['CONFIG']
assert sorted(tuple(r['patch']) for r in result['patches'])==[(i,j) for i in range(16) for j in range(i+1)]
for r in result['patches']:
 assert r['passed'] and not r['failed_leaves']
 assert r['nodes_tested']==1 and r['accepted_leaves'][0]['path']==''
 assert C['unpack'](r['accepted_leaves'][0]['H_lower'])>C['TARGET']
record={'status':'PASS finite exact/Arb selfchecks','exact_Bernstein_basis_checks':count,'DCT_identity_entries':N*N,'off_node_checks':rows,'all_patch_coverage_exact':True,'all_136_root_boxes_above_target':True,'seconds':time.perf_counter()-start,'checker_sha256':hashlib.sha256((OUT/'cc11_cheb_check.py').read_bytes()).hexdigest()}
(OUT/'cc11_cheb_selfcheck_result.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
