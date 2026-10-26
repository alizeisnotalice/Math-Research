"""Finite2048-parent zero-H independent integer-root budget consumer."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import gzip,json,hashlib,time,os
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_adaptive103');N=1<<160;BITS=384

def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def atomic(p,x):
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(serial(x)+'\n');os.replace(tmp,p)
def fifth_floor(n):
 if n==0:return 0
 x=1<<((n.bit_length()+4)//5)
 while True:
  y=(4*x+n//x**4)//5
  if y>=x:break
  x=y
 while (x+1)**5<=n:x+=1
 while x**5>n:x-=1
 assert x**5<=n<(x+1)**5;return x

def upperroot(q,k):
 den=1<<BITS;n=q.numerator*den**k//q.denominator;v=fifth_floor(n) if k==5 else isqrt(isqrt(n));a=F(v,den);assert a**k<=q;u=a if a**k==q else F(v+1,den);assert q<=u**k;return u

def budget(up,E,g):
 uh=F(up,N);B=upperroot(uh**4,5)+E;I=upperroot(B**5,4);C=upperroot((5*B/g)**5,4);margin=32**4*g**5-5**5*B**5;assert (B-E)**5>=uh**4 and I**4>=B**5 and C**4>=(5*B/g)**5;return {'UH':str(uh),'B_upper':str(B),'I_upper':str(I),'cost_upper':str(C),'cost_upper_display':float(C),'strict32':margin>0,'exact_32_gate_margin':str(margin),'valid_g':str(g),'no_tube_added':True,'candidate_pending_I_main_and_independent_qualification':True,'formal_C5_changed':False}

def snapshot(reg,E,g,event_rows):
 p=H/'progress.json'
 if not p.exists():return None
 first=p.read_bytes();progress=json.loads(first);statebytes=(S/'state.json.gz').read_bytes()
 if hashlib.sha256(statebytes).hexdigest()!=progress['state_sha256'] or first!=p.read_bytes():return None
 st=json.loads(gzip.decompress(statebytes));meta=st['adaptive103'];seq=meta['event_sequence'];assert seq==len(meta['parents_complete']);assert seq<=len(event_rows);eventhash='0'*64 if seq==0 else event_rows[seq-1]['event_hash'];assert eventhash==meta['event_hash'];nodes=st['nodes'];up=sum(nodes[k]['effU'] for k in st['roots']);assert up==meta['globalU_units'];b=budget(up,E,g)
 if not b['strict32']:return None
 ml=max(nodes[k]['level'] for k in st['frontier']);assert sum(nodes[k]['weight']*2**(4*(ml-nodes[k]['level'])) for k in st['frontier'])==2**(4*ml);out=S/'g105_first_strict32_state.json.gz';assert not out.exists();tmp=out.with_suffix('.gz.tmp');tmp.write_bytes(statebytes);os.replace(tmp,out);proof={'status':'numerical_g105_gate_positive_pending_I_independent_and_main','parent_count':seq,'event_hash':eventhash,'state_path':str(out),'state_sha256':sha(out),'source_g105_acceptance_sha256':reg['g105_acceptance_sha256'],'production_registration_sha256':reg['production_registration_sha256'],'budget':b,'weighted_volume':'1','physical_frontier':len(st['frontier']),'new_H_calls':0};atomic(H/'g105_first_strict32.json',proof);return proof

def main():
 began=time.monotonic();reg=json.loads((H/'budget_g105_registration.json').read_text());assert sha(__file__)==reg['consumer_sha256'] and sha(H/'registration.json')==reg['production_registration_sha256'];gp=Path(reg['g105_acceptance_path']);assert sha(gp)==reg['g105_acceptance_sha256'];ga=json.loads(gp.read_text());assert ga['status']=='primary_accepts_complete_same_model_global_g105' and F(ga['qualified_g'])==F(21,200000);prod=json.loads((H/'registration.json').read_text());E=F(prod['E076_upper']);g=F(21,200000);rows=[];prev='0'*64;first_positive=None;previousU=None;stable_snapshot=False
 while True:
  events=[];jp=S/'parent_events.jsonl'
  if jp.exists():
   for raw in jp.read_text().splitlines():
    try:events.append(json.loads(raw))
    except json.JSONDecodeError:break
  for ev in events[len(rows):]:
   d={k:v for k,v in ev.items() if k!='event_hash'};assert ev['seq']==len(rows)+1 and ev['previous_hash']==prev and hashlib.sha256(serial(d).encode()).hexdigest()==ev['event_hash'];up=ev['global_U_units'];assert previousU is None or up<=previousU;b=budget(up,E,g);rows.append({'parent_count':ev['seq'],'parent':ev['parent'],'event_hash':ev['event_hash'],**b});prev=ev['event_hash'];previousU=up
   if b['strict32'] and first_positive is None:first_positive=ev['seq']
  if first_positive is not None and not stable_snapshot:
   snap=snapshot(reg,E,g,events)
   if snap is not None:stable_snapshot=True;print('G105_FIRST_STRICT32_SNAPSHOT',serial(snap),flush=True)
  rc=json.loads((H/'run_receipt.json').read_text());done=rc['status']!='running'
  if done and rc['status']=='complete' and len(rows)<2048:continue
  if rows:
   atomic(H/'g105_budget_progress.json',{'status':'complete' if done else 'watching_finite2048_parent_checkpoints','pid':os.getpid(),'complete_parent_prefixes_checked':len(rows),'maximum_registered_parents':2048,'first_exact_positive_event':first_positive,'strict32_snapshot_saved':stable_snapshot,'latest':rows[-1],'new_H_calls':0})
  if done:
   assert len(rows)<=2048
   if rc['status']=='complete':assert len(rows)==2048 and F(rows[-1]['UH'])==F(rc['after']['UH_truncated_integral'])
   if first_positive is not None and not stable_snapshot:
    snap=snapshot(reg,E,g,events);assert snap is not None;stable_snapshot=True
   break
  time.sleep(20)
 out={'status':'complete' if rc['status']=='complete' else 'stopped_on_production_fault','production_status':rc['status'],'g105_acceptance_sha256':reg['g105_acceptance_sha256'],'g105':'21/200000','production_kept_g100':True,'all_complete_parent_prefixes_checked':len(rows),'first_exact_positive_event':first_positive,'strict32_snapshot_saved':stable_snapshot,'final':rows[-1] if rows else None,'all_roots_outward_by_exact_integer_powers':True,'new_H_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-began,'consumer_sha256':sha(__file__)};atomic(H/'g105_budget_result.json',out);tmp=H/'g105_prefix_budgets.json.gz.tmp';tmp.write_bytes(gzip.compress(serial(rows).encode(),mtime=0));os.replace(tmp,H/'g105_prefix_budgets.json.gz');print(serial(out),flush=True)
if __name__=='__main__':main()
