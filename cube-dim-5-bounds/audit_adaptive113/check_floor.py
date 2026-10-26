#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,gzip,os,time
P=Path(__file__).resolve().parent;T=Path('/Users/zhengzhihao/Desktop/T');H=T/'Math-Research-cube5-alt-cont/cube-dim-5-bounds/upper_adaptive103';N=1<<160

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while z:=f.read(1<<20):h.update(z)
 return h.hexdigest()
def read(p):
 p=Path(p);return json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else json.loads(p.read_text())
def save(p,x):
 p=Path(p);temp=Path(str(p)+'.tmp');temp.write_text(json.dumps(x,indent=2)+'\n');temp.replace(p)
def kth(q,k,bits=320):
 # Binary integer bracketing, independent of producer's Newton fifth-root routine.
 D=1<<bits;num=q.numerator*D**k;den=q.denominator;lo=0;hi=1<<((num.bit_length()+k-1)//k)
 while lo+1<hi:
  m=(lo+hi)//2
  if m**k*den<=num:lo=m
  else:hi=m
 a=F(lo,D);b=a if a**k==q else F(lo+1,D);assert a**k<=q<=b**k;return a,b

def main():
 t=time.monotonic();reg=read(P/'floor_registration.json');assert not (P/'floor_results.json').exists(),'already completed'
 assert sha(__file__)==reg['code_sha256']
 for r in reg['sources']:assert sha(r['path'])==r['sha256'],r['path']
 sc=read(H/'mean_floor_certificate_scope.json');st=read(sc['source_state_path']);assert sha(sc['source_state_path'])==sc['source_state_sha256'];nodes=st['nodes'];assert len(nodes)==97972 and len(st['frontier'])==94580
 means={}
 for r in sc['input_certificates']:
  assert sha(r['path'])==r['sha256']
  for c in read(r['path']):
   mi=c['mean_i'] if 'mean_i' in c else c['bounds']['M_i'];means[c['key']]=mi
 assert set(means)==set(nodes);own={}
 for index,(key,n) in enumerate(nodes.items(),1):
  mi=means[key];assert mi>0;a=isqrt(isqrt(mi**5));assert a**4<=mi**5<(a+1)**4;au=a+(a**4<mi**5);D=1<<(40+4*n['level']);w=n['weight'];lower=min(n['selfU'],w*a//D);upper=min(n['selfU'],-(-(w*au)//D));own[key]=(lower,upper)
  if index%8192==0:save(P/'floor_checkpoint.json',{'nodes_completed':index,'total':97972})
 memo={}
 def visit(key):
  if key in memo:return memo[key]
  lo,up=own[key];children=nodes[key]['children']
  if children:
   v=[visit(c) for c in children];lo=min(lo,sum(a for a,b in v));up=min(up,sum(b for a,b in v))
  assert lo<=up;memo[key]=(lo,up);return lo,up
 roots=[visit(k) for k in st['roots']];assert len(memo)==97972;zl,zu=F(sum(a for a,b in roots),N),F(sum(b for a,b in roots),N);actual=sc['exact_hypothetical_object'];assert zl==F(actual['hypothetical_UH_lower']) and zu==F(actual['hypothetical_UH_upper'])
 prod=read(H/'registration.json');E=F(prod['E076_upper']);Esource=read(T/'Math-Research-cube5-astra-upper/cube-dim-5-bounds/upper_shorttime076/results.json');assert prod['epsilon']=='1/512' and E==F(next(c for c in Esource['cases'] if c['epsilon']=='1/512')['E_norm']['upper'])
 gp=H/'sourceg105_primary_acceptance.json';assert sha(gp)=='b249cc86042f4f5a8040cfdf21dcb7eca132813d0e0d531d2a40db1627045744';g=F(read(gp)['qualified_g']);assert g==F(21,200000)
 BL=kth(zl**4,5)[0]+E;BU=kth(zu**4,5)[1]+E;IL=kth(BL**5,4)[0];IU=kth(BU**5,4)[1];TL,TU=kth(32**4*(g/5)**5,4)
 assert BL>E and (BL-E)**5<=zl**4 and IL**4<=BL**5
 assert (BU-E)**5>=zu**4 and IU**4>=BU**5
 assert IL>TU and TU-IL<0 and 5**5*BL**5>32**4*g**5
 assert F(actual['hypothetical_budget_I_lower'])<=IL and IU<=F(actual['hypothetical_budget_I_upper'])
 assert F(sc['threshold_T_g105_lower'])<=TL<=TU<=F(sc['threshold_T_g105_upper'])
 assert F(actual['hypothetical_S_upper'])==F(sc['threshold_T_g105_upper'])-F(actual['hypothetical_budget_I_lower'])<0
 gap=read(H/'gap_analysis.json');AL=F(gap['fullG_I_lower']);AU=F(gap['fullG_I_upper']);assert AL==F(sum(nodes[k]['effL'] for k in st['roots']),N)
 assert gap['certificate_S_positive']==False and F(gap['certificate_S_lower'])==F(sc['threshold_T_g105_lower'])-AU and F(gap['certificate_S_upper'])==F(sc['threshold_T_g105_upper'])-AU
 assert F(gap['fixed_g105_trueI_S_lower'])==F(sc['threshold_T_g105_lower'])-AU and F(gap['fixed_g105_trueI_S_upper'])==F(sc['threshold_T_g105_upper'])-AL>0
 cl,cu=kth((5*BL/g)**5,4)[0],kth((5*BU/g)**5,4)[1];assert cl>32
 result={'status':'pass','scope':'fixed-current-M_U same-fixed-E nonnegative-variance upper-certificate family only','not_actual_I_or_model_C5_lower_bound':True,'nodes':97972,'physical_leaves':94580,'hypothetical_UH_lower':str(zl),'hypothetical_UH_upper':str(zu),'independent_I_expression_lower':str(IL),'independent_I_expression_upper':str(IU),'T_lower':str(TL),'T_upper':str(TU),'family_S_lower':str(TL-IU),'family_S_upper':str(TU-IL),'family_cost_lower':str(cl),'family_cost_upper':str(cu),'strict_family_exclusion':True,'actual_certificate_S_negative':True,'true_I_fixed_g105_S_interval_crosses_zero':True,'new_H_Arb_log_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-t,'code_sha256':sha(__file__),'registration_sha256':sha(P/'floor_registration.json')};save(P/'floor_results.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
