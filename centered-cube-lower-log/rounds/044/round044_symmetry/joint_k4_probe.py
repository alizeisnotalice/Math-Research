"""Numerical true-source K4 coupled score/gate integral, not a bound certificate."""
from pathlib import Path
import importlib.util,json,math,time,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'round042_selector/probe.py'
if not SOURCE.exists():SOURCE=Path(__file__).resolve().parents[2]/'042/round042_selector/probe.py'
spec=importlib.util.spec_from_file_location('source42',SOURCE);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 outpath=HERE/'JOINT_K4.json';pf=HERE/'JOINT_K4_PLAN.json';assert not outpath.exists()and not pf.exists();start=time.time()
 cfg=dict(cases=[dict(n=1,grid=512),dict(n=2,grid=128),dict(n=3,grid=64)],lengths=[1.2,1.5,1.8],batch=2048,selector='sum(j*alpha_j) mod4')
 pf.write_text(json.dumps(dict(config=cfg,script_SHA256=sha(__file__),helper42_SHA256=sha(SOURCE)),indent=2)+'\n');rows=[]
 for case in cfg['cases']:
  n=case['n'];grid=case['grid'];total=grid**n;c=old.CONFIG.copy();c.update(n=n,lengths=cfg['lengths']);types=old.states(n,4)
  s=np.array([sum(j*a[j]for j in range(4))%4 for a in types]);pi=np.bincount(s,weights=old.weights(types,c),minlength=4);assert np.array_equal(pi,np.full(4,.25))
  sums=np.zeros((3,3));maxcs=-np.inf
  for at in range(0,total,cfg['batch']):
   ids=np.arange(at,min(at+cfg['batch'],total));x=np.empty((len(ids),n));v=ids.copy()
   for i in range(n):x[:,i]=(v%grid+.5)/grid;v//=grid
   feat,fixed,_=old.features(x,c);F=old.selector_value(feat,s);G=np.einsum('b,bnj->nj',pi,fixed);sep=np.zeros_like(F)
   for b in range(4):
    comp=feat[b][...,s==b].sum(-1);sep+=comp*comp/(pi[b]*fixed[b])
   maxcs=max(maxcs,float(np.max(F*F/G-sep)));sums+=np.array([(F*F/G).sum(0),sep.sum(0),F.sum(0)])
  rows.append(dict(n=n,grid=grid,total=total,whole_F2_over_G=(sums[0]/total).tolist(),CS_separated_integral=(sums[1]/total).tolist(),F_mass_quadrature=(sums[2]/total).tolist(),max_pointwise_CS_violation=maxcs))
 out=dict(config=cfg,cases=rows,plan_SHA256=sha(pf),seconds=time.time()-start,scope='Actual K4 common-type source from42; no outward quadrature certificate or Cn claim')
 outpath.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
