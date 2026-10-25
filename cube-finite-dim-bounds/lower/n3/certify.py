"""Fraction enclosure for one rationalized finite two-orbit witness union.
No floating point enters the verifier after reading exact decimal rational inputs.
"""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parent

def floor_root(n,k):
 lo,hi=0,1<<((n.bit_length()+k-1)//k)
 while lo+1<hi:
  m=(lo+hi)//2
  if m**k<=n:lo=m
  else:hi=m
 return lo

def root_interval(m,D=10**24):
 q=floor_root(m.numerator*D**3//m.denominator,3)
 assert q**3*m.denominator<=m.numerator*D**3<(q+1)**3*m.denominator
 return F(q,D)/2,F(q+1,D)/2

def volume(boxes):
 @lru_cache(None)
 def rec(B):
  if not B:return F(0)
  if len(B[0])==2:
   a,b=B[0];s=F(0)
   for c,e in B[1:]:
    if c>b:s+=b-a;a,b=c,e
    else:b=max(b,e)
   return s+b-a
  cuts=sorted({x for b in B for x in b[:2]})
  return sum(((v-u)*rec(tuple(sorted({b[2:] for b in B if b[0]<=u and b[1]>=v}))) for u,v in zip(cuts,cuts[1:])),F(0))
 return rec(tuple(sorted(set(boxes))))

def index_collar(inp):
 """Integer collar for translations of the full receiver cell [0,P].
 Orbit representatives retain their given phases, even outside [0,1).
 """
 active_blocks=[]
 for key in inp['witness_keys']:
  assert len(key)==12 and all(isinstance(k,int) for k in key)
  for orbit in range(2):
   lo=key[3*orbit:3*orbit+3];hi=key[6+3*orbit:9+3*orbit]
   if all(l<=h for l,h in zip(lo,hi)):active_blocks.append((lo,hi))
 assert active_blocks
 index_min=[min(lo[j] for lo,hi in active_blocks) for j in range(3)]
 index_max=[max(hi[j] for lo,hi in active_blocks) for j in range(3)]
 b=max(abs(x) for x in index_min+index_max)
 assert all(-b<=lo[j]<=hi[j]<=b for lo,hi in active_blocks for j in range(3))
 mass=sum(map(F,inp['weights']),F(0))
 return {'index_min':index_min,'index_max':index_max,'b':b,
  'collar_convention':'receive physical cells P*(q+[0,1]^3), q in [-N,N]^3; retain both original orbit representatives at every integer index k in [-N-b,N+b]^3',
  'mass_per_index_pair':str(mass),
  'finite_input_mass_formula':str(mass)+'*(2*N+2*b+1)^3',
  'receiver_union_volume_formula':'(2*N+1)^3*U_tau',
  'phase_reindexing':'none; negative phase coordinates are kept and atom pairs are grouped by the integer orbit index, not by physical-cell membership'}

def verify(inp):
 p=list(map(F,inp['periods']));ph=[list(map(F,q)) for q in inp['phases']];w=list(map(F,inp['weights']))
 assert len(p)==3 and len(ph)==len(w)==2 and all(x>0 for x in p+w)
 inner=[];outer=[];details=[]
 for k in inp['witness_keys']:
  lo=[k[:3],k[3:6]];hi=[k[6:9],k[9:12]]
  counts=[max(0,hi[b][0]-lo[b][0]+1)*max(0,hi[b][1]-lo[b][1]+1)*max(0,hi[b][2]-lo[b][2]+1) for b in range(2)]
  mass=sum((w[b]*counts[b] for b in range(2)),F(0));rlo,rhi=root_interval(mass)
  active=[b for b in range(2) if counts[b]>0]
  lc=[max((hi[b][j]+ph[b][j])*p[j] for b in active) for j in range(3)]
  hc=[min((lo[b][j]+ph[b][j])*p[j] for b in active) for j in range(3)]
  for target,r in [(inner,rlo),(outer,rhi)]:
   bb=[]
   for j in range(3):bb.extend([max(F(0),lc[j]-r),min(p[j],hc[j]+r)])
   if all(bb[2*j]<bb[2*j+1] for j in range(3)):target.append(tuple(bb))
  details.append({'mass':str(mass),'radius_interval':[str(rlo),str(rhi)],'strict_inner_radius':(2*rlo)**3<mass})
 vlow=volume(inner)/sum(w);vhigh=volume(outer)/sum(w)
 cert=json.loads((ROOT.parents[1]/'baseline/original_certificates.json').read_text())[0]
 baseline_lower,baseline_upper=map(F,cert['lower_bound_interval'])
 return {**index_collar(inp),'status':'complete_certification_of_finite_witness_union_only','ratio_interval':[str(vlow),str(vhigh)],'decimal_display_only':[float(vlow),float(vhigh)],'strictly_below_original_R3':vhigh<baseline_lower,'rigorous_gap_lower':str(baseline_lower-vhigh),'witness_count':len(inp['witness_keys']),'nonempty_inner_boxes':len(inner),'nonempty_outer_boxes':len(outer),'radius_enclosures':details,'critical_inner_radius_count':sum(not x['strict_inner_radius'] for x in details),'strict_superlevel_handling':'uniform inward contraction then slightly smaller radius; finite-cell truncation and same-mass L1 smoothing; take limits'}

def main():
 start=time.monotonic()
 ap=argparse.ArgumentParser();ap.add_argument('--rebuild-input',action='store_true');args=ap.parse_args()
 if args.rebuild_input:
  r=json.loads((ROOT/'refine_results.json').read_text());c=r['best']
  inp={'periods':[format(x,'.12f') for x in c['p']],'phases':[['0','0','0'],[format(x,'.12f') for x in c['phase']]],'weights':[format(c['weights'][0],'.12f'),str(F(2)-F(format(c['weights'][0],'.12f')))],'witness_keys':json.loads((ROOT/'witness_family.json').read_text())}
  (ROOT/'certificate_input.json').write_text(json.dumps(inp,indent=2)+'\n')
 inp=json.loads((ROOT/'certificate_input.json').read_text())
 result=verify(inp);result['elapsed_seconds']=time.monotonic()-start;result['input_sha256']=hashlib.sha256((ROOT/'certificate_input.json').read_bytes()).hexdigest()
 (ROOT/'certificate_result.json').write_text(json.dumps(result,indent=2)+'\n')
 print({k:v for k,v in result.items() if k!='radius_enclosures'})
if __name__=='__main__':main()
