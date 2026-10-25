"""Fixed 2440-cell generation and independent-from-LP exact certificate check."""
from fractions import Fraction as Q
from pathlib import Path
import argparse,json,hashlib,time
from flint import fmpq,arb,ctx as arbctx
import ci4b_affine as engine
TARGET=Q(19,100000)
REG_SHA='a5696336c394920573905c491ad38c789a8d499831cf6865b36cb5ccdca68d53'
CAP_ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds')
CAP_LOCK={'proofs/round011_face_caps.md':'b1564bf4c0cb4a6e2f72786ea5ffae30c922ad45c0cb9874aae6bfa55d4fa8de','audit/030_audit.md':'dcbc9f31561668e0b31922990b659786038f205cdb106230921ace53a72db640'}
def lock_caps():
 for p,h in CAP_LOCK.items():assert hashlib.sha256((CAP_ROOT/p).read_bytes()).hexdigest()==h,(p,'cap-source drift')
def cells():
 rr=engine.roots();result=[]
 assert [engine.volume(t) for t in rr]==[Q(1,36),Q(1,12),Q(1,18)]
 for t in engine.refine(rr[0],8):result.append((0,t))
 for case in (1,2):
  for k,ts in engine.shells(rr[case]):
   for t in ts:
    parentvol=engine.volume(t);children=engine.refine(t,max(0,8-2*k))
    assert sum(engine.volume(z) for z in children)==parentvol
    result.extend((case,z) for z in children)
 assert len(result)==2440<=4096
 assert sum(engine.volume(t) for c,t in result)==Q(1,6)-Q(5,36*2**84)
 return result

def cap():
 pi=arb.pi();D=1/(48*pi**4);B=arb(180224)/pi**3;e=arb(1)/2**27
 return arb(2).root(3)*4*pi*((D.root(3)**4)*e**3/3+arb(3)/5*(B.root(3)**4)*e.root(3)**5)/8

def make(path):
 start=time.time();arbctx.prec=192
 assert hashlib.sha256(Path(__file__).with_name('ci4b_registration.md').read_bytes()).hexdigest()==REG_SHA
 lock_caps();models=engine.models();cc=cells();rows=[];subtotal=arb(0)
 for i,(case,tet) in enumerate(cc):
  r=engine.majorant(models[case],tet)
  if r['status']!='exact_coefficient_majorant':raise RuntimeError((i,case,r))
  pp=[Q(x) for x in r['P']]
  val=engine.aq(engine.volume(tet))*sum(engine.aq(p).root(3)**4 if p else arb(0) for p in pp)/(4*arb.pi()**4)
  subtotal+=val
  rows.append({'case':case,'vertices':[[str(x) for x in v] for v in tet],'P':r['P']})
  if (i+1)%100==0:print('certified cells',i+1,'of',len(cc),'elapsed',round(time.time()-start,2),flush=True)
 budget=cap();total=6*subtotal+budget
 out={'status':'upper certificate candidate; exact checker required','registration_sha256':REG_SHA,'source_sha256':engine.SOURCE_SHA,'interface_sha256':engine.INTERFACE_SHA,'cell_count':len(cc),'ordered_covered_volume':str(Q(1,6)-Q(5,36*2**84)),'omitted_corner_radius':'1/'+str(2**27),'precision_bits':arbctx.prec,'target':str(TARGET),'I4_upper_expression_display':str(total),'I4_upper_rational':str(total.upper().fmpq()),'cap_upper_rational':str(budget.upper().fmpq()),'target_pass':Q(str(total.upper().fmpq()))<=TARGET,'cells':rows,'elapsed_seconds':time.time()-start}
 Path(path).write_text(json.dumps(out,separators=(',',':'))+'\n')
 print('FULL MAKE',out['target_pass'],out['I4_upper_expression_display'],'seconds',out['elapsed_seconds'],flush=True)

def check(path,save):
 start=time.time();arbctx.prec=256;raw=json.loads(Path(path).read_text())
 assert raw['registration_sha256']==REG_SHA
 assert raw['source_sha256']==engine.SOURCE_SHA and raw['interface_sha256']==engine.INTERFACE_SHA
 assert hashlib.sha256(Path(__file__).with_name('ci4b_registration.md').read_bytes()).hexdigest()==REG_SHA
 lock_caps();models=engine.models();cc=cells();assert len(raw['cells'])==len(cc)
 total=arb(0);minimum=None;constraints=0
 for i,((case,tet),row) in enumerate(zip(cc,raw['cells'])):
  assert case==row['case'] and [[str(x) for x in v] for v in tet]==row['vertices']
  pp=list(map(Q,row['P']));assert len(pp)==4 and all(p>=0 for p in pp)
  sub=[sum((engine.CTX.constant(fmpq(v[j].numerator,v[j].denominator))*x for v,x in zip(tet,engine.XS)),engine.CTX.constant(0)) for j in range(4)]
  N,D=models[case];N=N.compose(*sub);D=D.compose(*sub)
  assert D and all(c>=0 for c in D.to_dict().values())
  P=sum((engine.CTX.constant(fmpq(p.numerator,p.denominator))*x for p,x in zip(pp,engine.XS)),engine.CTX.constant(0))
  residual=P*D-N*sum(engine.XS,engine.CTX.constant(0));coeff=list(residual.to_dict().values())
  assert all(c>=0 for c in coeff),(i,'negative majorant coefficient')
  constraints+=len(coeff)
  if coeff:minimum=min(coeff) if minimum is None else min(minimum,min(coeff))
  total+=engine.aq(engine.volume(tet))*sum(engine.aq(p).root(3)**4 if p else arb(0) for p in pp)/(4*arb.pi()**4)
  if (i+1)%400==0:print('checked',i+1,'of',len(cc),'elapsed',round(time.time()-start,2),flush=True)
 budget=cap();total=6*total+budget;upper=Q(str(total.upper().fmpq()))
 result={'status':'EXACT CONTINUOUS I4 CERTIFICATE PASS' if upper<=TARGET else 'exact continuous upper bound misses target','source_sha256':engine.SOURCE_SHA,'interface_sha256':engine.INTERFACE_SHA,'certificate_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'cap_source_hashes':CAP_LOCK,'registration_sha256':REG_SHA,'cell_count':len(cc),'ordered_covered_volume':str(sum(engine.volume(t) for c,t in cc)),'omitted_ordered_volume':str(Q(5,36*2**84)),'I4_upper_rational':str(upper),'I4_upper_expression_display':str(total),'cap_upper_rational':str(budget.upper().fmpq()),'cap_display':str(budget),'target':str(TARGET),'exact_slack':str(TARGET-upper),'target_pass':upper<=TARGET,'precision_bits':arbctx.prec,'nonzero_nonnegative_majorant_coefficients':constraints,'min_nonzero_coefficient':str(minimum),'elapsed_seconds':time.time()-start}
 if save:Path(save).write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();group=p.add_mutually_exclusive_group(required=True);group.add_argument('--make');group.add_argument('--check');p.add_argument('--save');args=p.parse_args()
 if args.make:make(args.make)
 else:check(args.check,args.save)
