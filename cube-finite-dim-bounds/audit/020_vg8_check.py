"""Locked vg8 fixed inputs: independent parity/endpoint arithmetic only."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod
import json,hashlib,argparse
HERE=Path(__file__).resolve().parent
HASHES={'vg8_superperiod.md':'82e9a5628938bb50352b06dfad845ec1411ad62b3f759455c844b10aef44a857','vg8_parity.py':'42ac3828f335f75ec46c8c0ef7c6e5708895e304007b96284c95b79fb5a9ded8','vg8_blocks.json':'05b2f7804dfdf8d16749a844fe7ed0361379e1ac02ecfd1f714763f398ffb551'}
def atoms(lo,hi):return list(product(*(range(a,b+1) for a,b in zip(lo,hi))))
def split_formula(lo,hi):
    lows=[];highs=[]
    for parity in (0,1):
      a=-((-(lo[0]-parity))//2);b=(hi[0]-parity)//2
      if a>b or any(x>y for x,y in zip(lo,hi)):lows.append((0,)*4);highs.append((-1,)*4)
      else:lows.append((a,)+tuple(lo[1:]));highs.append((b,)+tuple(hi[1:]))
    return tuple(z for row in lows+highs for z in row)
def run(root):
    for name,h in HASHES.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h
    cfg=json.loads((root/'vg8_blocks.json').read_text());author=json.loads((root/'vg8_results.json').read_text());rows=[];total=0
    for block,a in zip(cfg['blocks'],author['rows']):
      lo,hi=block['lo'],block['hi'];raw=atoms(lo,hi);counts=[sum((k[0]%2)==p for k in raw) for p in (0,1)];N=len(raw);delta=counts[0]-counts[1];key=split_formula(lo,hi);constructed=[]
      for p in (0,1):
        for m in atoms(key[p*4:p*4+4],key[8+p*4:12+p*4]):constructed.append((2*m[0]+p,)+m[1:])
      assert set(constructed)==set(raw) and len(constructed)==len(raw)
      assert list(key)==a['key'] and counts==a['counts'] and N==a['N'] and delta==a['delta']
      T=prod(hi[j]-lo[j]+1 for j in range(1,4));assert delta==(0 if (hi[0]-lo[0]+1)%2==0 else (-1)**lo[0]*T)
      masses=[]
      for e in map(F,cfg['finite_mass_parameters']):
        exact=sum((1+e if k[0]%2==0 else 1-e for k in raw),F(0));assert exact==N+e*delta;masses.append(str(exact))
      assert masses==a['masses'];shiftlo=[lo[0]+1]+lo[1:];shifthi=[hi[0]+1]+hi[1:];shift=atoms(shiftlo,shifthi);assert [sum(k[0]%2==p for k in shift) for p in (0,1)]==counts[::-1]
      total+=N;rows.append({'counts':counts,'N':N,'delta':delta,'key':list(key)})
    assert len(rows)==12 and total==391
    eta=F(cfg['eta']);r0=F(3,2);q=F(4095,4096);tau=3*(1-q**4);rt=3*q/2;r1=F(-1,8);rem=rt-r0-r1*tau;c=F(-1,3);rho=abs(c)*eta;Br=3*r0*c*c/(32*(1-rho)**2)
    assert 0<tau<eta and (81-27*tau)/16==rt**4 and rem==F(-301940739,2251799813685248) and Br==F(576,36481) and abs(rem)<=Br*tau*tau
    # Reconstruct raw normalized source atoms from the first N81 block.
    base=cfg['blocks'][-2];key=split_formula(base['lo'],base['hi']);pts=[];basepts=[];vels=[]
    for p in (0,1):
      for m in atoms(key[p*4:p*4+4],key[8+p*4:12+p*4]):
        b=tuple(F(m[j])+(F(p,2) if j==0 else F(0)) for j in range(4));v=(F(0),F(p,4),F(0),F(0));basepts.append(b);vels.append(v);pts.append(tuple(b[j]+tau*v[j] for j in range(4)))
    Q=[2,1,1,1];errors=[];derivatives=[]
    for j in range(4):
      xmax=max(b[j] for b in basepts);xmin=min(b[j] for b in basepts);maxv=max(v[j] for b,v in zip(basepts,vels) if b[j]==xmax);minv=min(v[j] for b,v in zip(basepts,vels) if b[j]==xmin)
      for actual,b0,v0 in [(max(b[j] for b in pts)-rt/Q[j],xmax-r0/Q[j],maxv-r1/Q[j]),(min(b[j] for b in pts)+rt/Q[j],xmin+r0/Q[j],minv+r1/Q[j])]:
        error=abs(actual-b0-tau*v0);assert error<=Br*tau*tau/Q[j];errors.append(str(error));derivatives.append(str(v0))
    assert errors==author['finite_negative_mass_example']['endpoint_absolute_errors'] and derivatives[2:4]==['3/8','-1/8']
    return {'status':'PASS independent locked vg8 fixed interfaces','source_sha256':HASHES,'fixed_blocks':12,'enumerated_atom_occurrences':total,'rows':rows,'negative_slope_N81':{'tau':str(tau),'r0':str(r0),'r1':str(r1),'remainder':str(rem),'Br':str(Br),'endpoint_velocities':derivatives,'endpoint_absolute_errors':errors},'scope':'parity/mass/raw receiver/conditional periodic embedding interface only; no actual R4 union recomputation or new L4 certificate','shared':'Python int/Fraction; no author parity/guard/checker imports'}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,default=HERE.parent/'lower/variation_geometry');ap.add_argument('--save',action='store_true');args=ap.parse_args();out=run(args.source_dir);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:HERE.joinpath('020_vg8_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
