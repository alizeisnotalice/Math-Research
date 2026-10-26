#!/usr/bin/env python3
"""Full independent largest-coordinate moments, alternate R/mean/chord, hierarchy."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import gzip,json,hashlib,multiprocessing as mp,os,time,sys
from flint import arb,ctx,fmpq
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_reenvelope088');N=1<<160;TS=None

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def rq(x):
 x=F(str(x));return fmpq(x.numerator,x.denominator)
def ep(x,up=False):return F(str((x.upper() if up else x.lower()).fmpq()))
def multiply(p,q):
 out=[F(0)]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):out[i+j]+=a*b
 return out
def conditional_maximum_moment(spans):
 # Condition on the unique continuous coordinate attaining largest absolute value.
 # Distinct path from085's survival-CDF integral; all present widths are positive.
 total=F(0);knots={F(0)}|{z for a,b in spans for v in (a,b) for z in (abs(v),-abs(v))}
 for i,(lo,hi) in enumerate(spans):
  assert lo<hi;cuts=sorted({lo,hi}|{x for x in knots if lo<x<hi})
  for left,right in zip(cuts,cuts[1:]):
   mid=(left+right)/2;sgn=1 if mid>=0 else -1;r=abs(mid);poly=[F(1)]
   for j,(a,b) in enumerate(spans):
    if i==j:continue
    top=(b,F(0)) if b<r else (F(0),F(sgn));bottom=(a,F(0)) if a>-r else (F(0),F(-sgn));c0=top[0]-bottom[0];c1=top[1]-bottom[1]
    if c0+c1*mid<=0:poly=[F(0)];break
    poly=multiply(poly,[c0/(b-a),c1/(b-a)])
   contribution=sum((c*(right**(k+3)-left**(k+3))/((hi-lo)*(k+3)) for k,c in enumerate(poly)),F(0));assert contribution>=0;total+=contribution
 assert 0<=total<=max(abs(x) for span in spans for x in span)**2
 return total
def init(ts):
 global TS
 TS=ts;ctx.prec=384
def verify_one(row):
 rank=row['rank'];op=S/'nodes'/f'{rank//1024:02d}'/f'{rank:05d}.json.gz';r=json.loads(gzip.decompress(op.read_bytes()));assert r['status']=='complete' and r['key']==row['key'] and r['source_evidence_sha256']==row['source_evidence_sha256'];assert sha(row['source_evidence_path'])==row['source_evidence_sha256'];old=json.loads(gzip.decompress(Path(row['source_evidence_path']).read_bytes()));p=list(map(F,old['point'] if 'point' in old else old['anchor']));L=row['level'];spans=[(F(i,2**L)-x,F(i+1,2**L)-x) for i,x in zip(row['indices'],p[1:])];M2=conditional_maximum_moment(spans);assert M2==F(r['M2']);mu=[(a+b)/2 for a,b in spans];eta=F(2,9)*sum((b-a)**2/12 for a,b in spans)+sum(x*x for x in mu)/4-sum(mu)**2/36;hinf=max(abs(x) for span in spans for x in span);B2=max(sum(x*x for x in u)/4-sum(u)**2/36 for u in product(*spans));assert eta==F(r['eta_old_Q']) and hinf==F(r['h_infinity']) and B2==F(old['geometry']['B2']);lam=max(map(abs,mu));ends=old['endpoints'];assert [x['t'] for x in ends]==TS;mean=arb(0)
 for i,e in enumerate(ends):
  a=F(e['t']);nxt=ends[i+1]['H'] if i+1<73 else {'lower':'0','upper':'0'};ml=F(e['H']['lower'])-F(nxt['upper']);mh=F(e['H']['upper'])-F(nxt['lower']);assert 0<ml<=mh;z=arb(rq(1+B2/(10*a*a)));P=16*z**5-20*z**3+5*z;oldR=P+(P*P-1).sqrt();ex=arb(rq(hinf/(2*a))).exp()**2;R=min(ep(oldR,True),ep(ex,True));C=min(R,1+lam/a+F(25,48)*M2*R/a**2,1+lam/a+F(21,8)*eta*R/a**2);mean+=arb(rq(mh*C))
 b=old['bounds'];l=F(b['L_i'],N);u=F(b['U_i'],N);m=F(r['new_mean_i'],N);assert min(F(b['M_i'],N),ep(mean,True))<=m<=F(b['M_i'],N) and l<=m<=u;v=arb(1)/2**(4*L);q=arb(5)/4;lp=arb(rq(l))**q;up=arb(rq(u))**q;chord=v*(lp*arb(rq((u-m)/(u-l)))+up*arb(rq((m-l)/(u-l))));assert ep(chord,True)<=F(r['new_chord_integral_upper_i'],N);assert r['new_self_U']==min(row['old_self_U'],row['weight']*r['new_chord_integral_upper_i']) and row['old_self_L']<=r['new_self_U'];assert r['new_H_calls']==0 and r['positive_bands']==73
 return {'rank':rank,'key':row['key'],'record_sha256':hashlib.sha256(serial(r).encode()).hexdigest(),'pid':os.getpid()}
def hierarchy_readback(proofs,reg):
 initial=json.loads(gzip.decompress((S/'initial_state.json.gz').read_bytes()));state=initial;nodes=state['nodes'];front=list(state['frontier']);LH0=sum(nodes[k]['effL'] for k in state['roots']);proofmap={p['rank']:p for p in proofs};allrows=[];previous=None;previousU=F(reg['initial_UH'])
 for bid in range(14):
  tag=f'{bid:04d}';brp=H/('block'+tag+'_registration.json');br=json.loads(brp.read_text());rp=H/('block'+tag+'_receipt.json');rc=json.loads(rp.read_text());cp=H/('block'+tag+'_cells.json.gz');assert sha(cp)==rc['certificates_sha256'];assert sha(brp)==rc['registration_sha256'];assert br['parent_receipt_sha256']==previous;rows=json.loads(gzip.decompress(cp.read_bytes()));assert len(rows)==br['count']<=4096;assert [r['rank'] for r in rows]==list(range(br['start'],br['end']));assert [r['key'] for r in rows]==br['keys']
  for r in rows:
   assert hashlib.sha256(serial(r).encode()).hexdigest()==proofmap[r['rank']]['record_sha256'];n=nodes[r['key']];assert n['selfL']==r['old_self_L'] and n['selfU']==r['old_self_U'];n['selfU']=r['new_self_U']
  def visit(k):
   n=nodes[k]
   if n['children']:
    for c in n['children']:visit(c)
    n['effU']=min(n['selfU'],sum(nodes[c]['effU'] for c in n['children']));n['effL']=max(n['selfL'],sum(nodes[c]['effL'] for c in n['children']))
   else:n['effU'],n['effL']=n['selfU'],n['selfL']
   assert n['effL']<=n['effU']
  for k in state['roots']:visit(k)
  lo=F(sum(nodes[k]['effL'] for k in state['roots']),N);up=F(sum(nodes[k]['effU'] for k in state['roots']),N);assert lo==F(reg['initial_LH']) and sum(nodes[k]['effL'] for k in state['roots'])==LH0;assert up<=previousU and up==F(rc['after']['UH_truncated_integral']);E=F(reg['E076_upper']);B=ep(arb(rq(up))**(arb(4)/5)+arb(rq(E)),True);assert B<=F(rc['after']['B_upper']);g=F(reg['valid_g']);assert (5**5*F(rc['after']['B_upper'])**5<32**4*g**5)==rc['after']['cost_strictly_below32'];assert state['frontier']==front;previousU=up;previous=sha(rp);allrows+=rows
 maxlevel=max(nodes[k]['level'] for k in front);assert sum(nodes[k]['weight']*2**(4*(maxlevel-nodes[k]['level'])) for k in front)==2**(4*maxlevel);actual=json.loads(gzip.decompress((S/'state.json.gz').read_bytes()));assert actual['nodes']==nodes and actual['frontier']==front;return up,allrows

def main():
 t=time.monotonic();ctx.prec=384;reg=json.loads((H/'registration.json').read_text());assert sha(H/'engine.py')==reg['engine_sha256'] and sha(H/'moment085.py')==reg['moment085_sha256'];assert sha(reg['g_acceptance_path'])==reg['g_acceptance_sha256'];ga=json.loads(Path(reg['g_acceptance_path']).read_text());assert ga['g']=='1/10000' and ga['status']=='primary_accepts_complete_alpha1_cube_g_component';rows=json.loads(gzip.decompress((H/'nodes.json.gz').read_bytes()));assert len(rows)==56139;pool=mp.get_context('spawn').Pool(12,initializer=init,initargs=(reg['time_endpoints'],));proofs=[]
 try:
  for p in pool.imap_unordered(verify_one,rows):
   proofs.append(p)
   if len(proofs)%8192==0:print('VERIFIED',len(proofs),round(time.monotonic()-t,3),flush=True)
  pool.close();pool.join()
 except BaseException:pool.terminate();pool.join();raise
 up,certs=hierarchy_readback(proofs,reg);result={'status':'pass','nodes_verified':len(proofs),'all_conditional_maximum_moments_match085':True,'all73_saved_band_masses_positive':True,'positive_bands_checked':73*len(proofs),'alternate_R_mean_chord_checked':True,'self_lowers_unchanged':True,'all14_hierarchy_prefixes_valid':True,'physical_frontier_unchanged':55819,'weighted_volume':'1','valid_g':reg['valid_g'],'g_acceptance_sha256':reg['g_acceptance_sha256'],'UH':str(up),'new_H_calls':0,'all12_workers_joined':True,'worker_pids':sorted({p['pid'] for p in proofs}),'pid':os.getpid(),'seconds':time.monotonic()-t,'source_sha256':sha(__file__)};(H/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
