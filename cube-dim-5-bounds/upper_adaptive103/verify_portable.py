"""Portable integer-only ledger reconstruction; no raw H or external scratch needed."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product
from math import factorial,prod
import gzip,json,hashlib,time,os
H=Path(__file__).resolve().parent;B=H.parent;N=1<<160

def read(p):return json.loads(gzip.decompress(p.read_bytes())) if str(p).endswith('.gz') else json.loads(p.read_text())
def key(L,i):return str(L)+':'+','.join(map(str,i))
def weight(i):return factorial(4)//prod(factorial(v) for v in Counter(i).values())
def serial(x):return json.dumps(x,sort_keys=True,separators=(',',':'))
def rebuild095():
 roots=read(B/'upper_global080/roots.json.gz');table=read(B/'upper_global080/initial_cells.json.gz');nodes={}
 for xs in table['rows']:
  r=dict(zip(table['columns'],xs));ii=roots[r['rank']]['indices'];k=key(5,ii);nodes[k]={'key':k,'level':5,'indices':ii,'weight':r['weight'],'selfL':r['weight']*r['integral_lower_i'],'selfU':r['weight']*r['integral_upper_i'],'children':None,'parent':None}
 rootkeys=sorted(nodes);front=set(nodes)
 for bid in range(2):
  tag=f'{bid:04d}';reg=read(B/f'upper_adapt082/batch{tag}_registration.json');cert={r['key']:r for r in read(B/f'upper_adapt082/batch{tag}_cells.json.gz')}
  for p,cs in reg['children_by_parent'].items():
   assert p in front;par=nodes[p];assert set(cs)=={key(par['level']+1,sorted((2*i+b for i,b in zip(par['indices'],bits)),reverse=True)) for bits in product((0,1),repeat=4)}
   for c in cs:
    r=cert[c];bb=r['bounds'];nodes[c]={'key':c,'level':r['level'],'indices':r['indices'],'weight':r['weight'],'selfL':r['weight']*bb['integral_lower_i'],'selfU':r['weight']*bb['integral_upper_i'],'children':None,'parent':p}
   par['children']=cs;front.remove(p);front.update(cs)
 for bid in range(14):
  for r in read(B/f'upper_reenvelope088/block{bid:04d}_cells.json.gz'):
   n=nodes[r['key']];assert n['selfL']==r['old_self_L'] and n['selfU']==r['old_self_U'];n['selfU']=r['new_self_U']
 for bid in range(14):
  for r in read(B/f'upper_range092/block{bid:04d}_cells.json.gz'):
   n=nodes[r['key']];assert n['selfL']==r['old_self_L'] and n['selfU']==r['old_self_U'];n['selfL']=r['new_self_L'];n['selfU']=r['new_self_U']
 def visit(k):
  n=nodes[k]
  if n['children']:
   for c in n['children']:visit(c)
   n['effL']=max(n['selfL'],sum(nodes[c]['effL'] for c in n['children']));n['effU']=min(n['selfU'],sum(nodes[c]['effU'] for c in n['children']))
  else:n['effL'],n['effU']=n['selfL'],n['selfU']
  assert n['effL']<=n['effU']
 for k in rootkeys:visit(k)
 # Continue through the already frozen095 complete-parent certificates.
 p095=B/'upper_adaptive095';r095=read(p095/'preflight_registration.json');c095={r['key']:r for r in read(p095/'children_certificates.json.gz')}
 for p,cs in r095['children_by_parent'].items():
  assert p in front;par=nodes[p]
  for c in cs:
   r=c095[c];bb=r['bounds'];nodes[c]={'key':c,'level':r['level'],'indices':r['indices'],'weight':r['weight'],'selfL':r['weight']*bb['integral_lower_i'],'selfU':r['weight']*bb['integral_upper_i'],'children':None,'parent':p}
  assert sum(nodes[c]['weight'] for c in cs)==16*par['weight'];par['children']=cs;front.remove(p);front.update(cs)
 for k in rootkeys:visit(k)
 a095=read(p095/'run_receipt.json')['after'];assert F(sum(nodes[k]['effL'] for k in rootkeys),N)==F(a095['LH_truncated_integral']) and F(sum(nodes[k]['effU'] for k in rootkeys),N)==F(a095['UH_truncated_integral'])
 return nodes,rootkeys,front

def main():
 t=time.monotonic();nodes,roots,front=rebuild095();pre=read(H/'preflight_registration.json');expected=sorted(front,key=lambda k:(-(nodes[k]['selfU']-nodes[k]['selfL']),nodes[k]['level'],nodes[k]['indices']))[:2048];assert expected==[x['key'] for x in pre['parents']];cert={c['key']:c for c in read(H/'children_certificates.json.gz')};jobs={j['key']:j for j in read(H/'jobs.json.gz')};receipt=read(H/'run_receipt.json');assert sum(nodes[k]['effL'] for k in roots)==F(receipt['before']['LH_truncated_integral'])*N;assert sum(nodes[k]['effU'] for k in roots)==F(receipt['before']['UH_truncated_integral'])*N
 for p in expected:
  par=nodes[p];assert p in front;cs=pre['children_by_parent'][p];assert set(cs)=={key(par['level']+1,sorted((2*i+b for i,b in zip(par['indices'],bits)),reverse=True)) for bits in product((0,1),repeat=4)};assert sum(cert[c]['weight'] for c in cs)==16*par['weight']
  for c in cs:
   r=cert[c];j=jobs[c];b=r['bounds'];ii=r['indices'];L=r['level'];assert r['weight']==weight(ii) and [i//2 for i in ii]==par['indices'] and j['parent']==p;l,u,m=F(b['L_i'],N),F(b['U_i'],N),F(b['M_i'],N);assert 0<l<=m<=u<=3;v=F(1,2**(4*L));assert (F(b['integral_lower_i'],N)/v)**4<=l**5
   # Convex-chord comparison using exact algebra and an independently rounded root.
   def root(q):
    from math import isqrt
    den=1<<224;n=isqrt(isqrt(q.numerator*den**4//q.denominator));return F(n+1,den)
   cl=v*root(l**5) if l==u else v*((u-m)*root(l**5)+(m-l)*root(u**5))/(u-l);assert cl<=F(b['integral_upper_i'],N);nodes[c]={'key':c,'level':L,'indices':ii,'weight':r['weight'],'selfL':r['weight']*b['integral_lower_i'],'selfU':r['weight']*b['integral_upper_i'],'effL':r['weight']*b['integral_lower_i'],'effU':r['weight']*b['integral_upper_i'],'children':None,'parent':p}
  par['children']=cs;front.remove(p);front.update(cs)
 for k in sorted(nodes,key=lambda k:nodes[k]['level'],reverse=True):
  n=nodes[k]
  if n['children']:n['effL']=max(n['selfL'],sum(nodes[c]['effL'] for c in n['children']));n['effU']=min(n['selfU'],sum(nodes[c]['effU'] for c in n['children']))
  assert n['effL']<=n['effU']
 lo,up=F(sum(nodes[k]['effL'] for k in roots),N),F(sum(nodes[k]['effU'] for k in roots),N);b=receipt['after'];assert lo==F(b['LH_truncated_integral']) and up==F(b['UH_truncated_integral']);assert sorted(front)==read(H/'physical_frontier.json.gz');ml=max(nodes[k]['level'] for k in front);assert sum(nodes[k]['weight']*2**(4*(ml-nodes[k]['level'])) for k in front)==2**(4*ml);reg=read(H/'registration.json');Bv,E,I,C,g=map(F,[b['B_upper'],reg['E076_upper'],b['I_upper'],b['cost_upper'],reg['valid_g']]);assert (Bv-E)**5>=up**4 and I**4>=Bv**5 and C**4>=(5*Bv/g)**5;assert (32**4*g**5>5**5*Bv**5)==b['strict32'];out={'status':'pass','all_source080_082_088_092_095_ledgers_reconstructed':True,'children':len(cert),'nodes':len(nodes),'physical_frontier':len(front),'weighted_volume':'1','LH':str(lo),'UH':str(up),'strict32':b['strict32'],'new_H_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-t,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(H/'portable_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
