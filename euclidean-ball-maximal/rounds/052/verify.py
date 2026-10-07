#!/usr/bin/env python3
"""Recompute the actual-source profiles, cap kernels and cutoff consistency checks."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
import hashlib,json,os,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[1]
def cutoff_checks():
 rows=[];checks=0
 for batch,ns in enumerate(((1,2,8,16,32),(64,128,256,512,1024),(2048,4096,8192,16384,32768))):
  for n in ns:
   B=n*(2*n-1).bit_length();L=0
   while 2**(L+1)<B+2:L+=1
   if n>=32:
    assert 7<=L<=n and F(4,3*L)*(1+F(64,3*L))<=F(340,441)<F(4,5)
    assert F(4,5)**4<F(1,2)
    assert F(B,2)+1<=2**L
   for s in (0,1,2,4,8):
    if n<32:
     bound=2*(n-1)*F(2)**(-n-s);assert bound<=F(2)**(-s-1)
    else:
     with localcontext() as ctx:
      ctx.prec=80;log2=D(2).ln();rho=(-D(L+s)*log2/n).exp();delta=1-rho*rho
      K=1+D(n)/2*((2/delta).ln()+(1+16*rho*rho/(n*delta)).ln())/log2
      bound=(D(2)**(-L-s))*(K+1)
      assert K<=D(B)/2 and bound<=D(2)**(-s)
    checks+=1
   rows.append(dict(batch=batch,n=n,B=B,base_depth=L if n>=32 else n,base_radius='2^(-L/n)' if n>=32 else '1/2',tested_extra_depths=[0,1,2,4,8]))
 return dict(status='passed',checks=checks,batch_counts=[5,5,5],cases=rows,scope='Exact integer constant chain; 80-digit Decimal envelope checks are not interval certificates. General summability follows analytically.')
def main():
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 with tempfile.TemporaryDirectory(prefix='euclidean52_') as folder:
  data={}
  for name in ('source','geometry'):
   out=Path(folder)/(name+'.json')
   subprocess.run([sys.executable,str(HERE/(name+'_probe.py')),'--output',str(out)],check=True,capture_output=True,text=True,env=env)
   data[name]=json.loads(out.read_text())
  source=data['source'];geo=data['geometry']
  assert source['status']==geo['status']=='passed' and source['batch_counts']==[3,3,3]
  assert len(source['cases'])==9 and sum(c['source_terms'] for c in source['cases'])==1463
  assert all(c['independent_round49_O_check']=='passed' for c in source['cases'])
  assert geo['geometry_case_count']==25
  assert sum(r['exact_shell_checks'] for r in geo['rational_kernel'])==1053
  assert sum(r['identity_checks'] for r in geo['rational_kernel'])==27
  keys=['batch','n','kind','parameter','terms','A','G','K','cross_precision_max_absolute_difference']
  compact=[{key:c[key] for key in keys}|{key:c[key] for key in ['B','L','rho_n_K_plus_one','lower_bound','exact_G'] if key in c} for c in geo['geometry']]
  deps=sorted((PROJECT/'runtime').rglob('*.py'))+[PROJECT/('rounds/'+p) for p in ['041/square_probe.py','042/row_probe.py','043/global_probe.py','044/tail_probe.py','046/prefix_probe.py','046/threshold_probe.py','047/flow_probe.py','049/global_probe.py','049/verification.json']]+sorted(HERE.glob('*.py'))
  receipt=dict(status='passed',round=52,date='2026-10-07',python=sys.version.split()[0],source=source,geometry=dict(status='passed',dimensions=geo['dimensions'],case_count=25,rational_kernel=geo['rational_kernel'],cases=compact,scope='120/160-digit comparison, not certified error bounds; masks not asserted to be actual exits.'),cutoff=cutoff_checks(),main_theorem_proved=False,boundary_shell_budget_proved=False,sha256={str(p.relative_to(PROJECT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in deps})
  print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
