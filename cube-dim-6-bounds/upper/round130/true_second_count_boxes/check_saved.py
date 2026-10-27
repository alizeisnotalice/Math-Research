#!/usr/bin/env python3
"""Saved leaves only: rational coverage/sums and fresh Arb768 corner enclosures."""
import sys
sys.set_int_max_str_digits(0)
import gzip,hashlib,json,math
from pathlib import Path
from fractions import Fraction as Q
from flint import arb,fmpq,ctx
ctx.prec=768
R=Path(__file__).resolve().parent
def rd(n):return json.loads((R/n).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def A(q):q=Q(q);return arb(fmpq(q.numerator,q.denominator))
def ends(value):
    assert value.is_finite();out=[]
    for z in (value.lower(),value.upper()):m,e=z.man_exp();out.append(Q(int(m))*Q(2)**int(e))
    return out
def pairs(v):return sum((v[i]*v[j]for i in range(32)for j in range(i+1,32)),Q(0))
spec=rd('inputs.json');table=rd('cone_table.json');result=rd('results.json');assert result['status']=='COMPLETE_FIXED_FINITE_BOX_PILOT'
for s in rd('source_manifest.json')['sources']:assert sha(s['path'])==s['sha256']
pL,pU=map(Q,spec['p_interval']);h=Q(spec['h']);sh=A(h).sqrt();ts=1-Q(1,2**24);tailprob=1-ts**6
fl,fu=map(Q,result['pi_factor_interval']);pl,pu=ends(A(46080)/arb.pi()**3);assert fl<=pl<=pu<=fu
assert len(result['leaf_files'])==len(table['cones'])==21
assert sorted(c['id']for c in table['cones'])==list(range(21))
assert len({entry['path']for entry in result['leaf_files']})==21
low=up=boxwidth=Q(0);kernels=0;checked_boxes=0;leafcounts=[];boxcounts=[];expected_producer_kernel_calls=21*64;cone_ids=[];root_low=Q(0)
for cid,entry in enumerate(result['leaf_files']):
 p=R/entry['path'];assert sha(p)==entry['sha256']
 with gzip.open(p,'rt')as f:data=json.load(f)
 cone=table['cones'][cid];assert data['cone_id']==cid==cone['id'] and data['cone_table_sha256']==sha(R/'cone_table.json')
 rows=data['leaves'];assert 1<=len(rows)<=256 and data['evaluated_boxes']==2*len(rows)-1
 cone_ids.append(data['cone_id']);leafcounts.append(len(rows));boxcounts.append(2*len(rows)-1);vertices=set();root_upper_corner=None
 prefixes=sorted(row['prefix']for row in rows)
 assert len(prefixes)==len(set(prefixes)) and sum((Q(1,2**len(z))for z in prefixes),Q(0))==1
 assert not any(b.startswith(a)for a,b in zip(prefixes,prefixes[1:]))
 volume=cl=cu=Q(0);nodes=[[Q(v)for v in row]for row in cone['nodes_by_body']]
 for row in rows:
  lo=[Q(0)]*6;hi=[ts]*6
  for digit in row['prefix']:
   assert digit in '01';axis=max(range(6),key=lambda j:(hi[j]-lo[j],-j));mid=(lo[axis]+hi[axis])/2
   if digit=='0':hi[axis]=mid
   else:lo[axis]=mid
  assert lo==list(map(Q,row['lower_t']))and hi==list(map(Q,row['upper_t']))
  vertices.update((tuple(lo),tuple(hi)))
  if hi==[ts]*6:root_upper_corner=row['kernel_at_upper_t']
  vol=math.prod(b-a for a,b in zip(lo,hi));assert vol==Q(row['volume']);volume+=vol;corner=[]
  for t,saved in [(lo,row['kernel_at_lower_t']),(hi,row['kernel_at_upper_t'])]:
   E=[-(1-A(v)).log()for v in t];iv=[];assert len(saved)==32
   for b,pair in zip(nodes,saved):
    if all(z==0 or t[j]==0 for j,z in enumerate(b)):value=sh.erf()
    else:
     a=sum((A(z)*E[j]for j,z in enumerate(b)),arb(0));u=a/(2*sh);value=((u-sh).erfc()-(2*a).exp()*(u+sh).erfc())/2
    l,u=ends(value);sl,su=map(Q,pair);assert 0<=sl<=max(Q(0),l)<=min(pU,u)<=su<=pU
    iv.append([sl,su]);kernels+=1
   corner.append(iv)
  lower=pairs([v[0]for v in corner[1]]);upper=pairs([v[1]for v in corner[0]])
  assert list(map(Q,row['F_interval']))==[lower,upper]
  assert list(map(Q,row['integral_interval']))==[vol*lower,vol*upper]
  cl+=vol*lower;cu+=vol*upper;checked_boxes+=1
 assert volume==ts**6 and root_upper_corner is not None
 # Every discarded ancestor corner survives as a final-leaf corner. With at most
 #512 vertices/cone there is no eviction from the2048-entry cache; bind clears it.
 assert len(vertices)<=2*len(rows)<=512
 if len(rows)>1:expected_producer_kernel_calls+=32*len(vertices)
 root_low+=Q(cone['alpha'])*ts**6*pairs([Q(v[0])for v in root_upper_corner])
 boxwidth+=Q(cone['alpha'])*(cu-cl)
 low+=Q(cone['alpha'])*cl;up+=Q(cone['alpha'])*(cu+496*pU**2*tailprob)
assert sorted(cone_ids)==list(range(21)) and len(set(cone_ids))==21
assert result['leaf_counts']==leafcounts and result['evaluated_boxes_by_cone']==boxcounts
assert result['evaluated_boxes']==sum(boxcounts)<=10731
assert result['kernel_cache_per_cone_maxsize']==2048
assert result['kernel_calls']==expected_producer_kernel_calls<=686784
assert result['precision_bits']==512 and result['h']=='1/128'and result['tail_L']==24
assert result['pair_terms']==496 and result['cone_count']==21
assert [46080*low,46080*up]==list(map(Q,result['pi_cubed_weighted_interval']))
assert [fl*low,fu*up]==list(map(Q,result['moment_interval']))
oldL,oldU=map(Q,spec['old_beta2_interval']);assert result['strict_saved128_beta2_improvement']==(fu*up<oldU)
assert result['strict_old_cap_expression_improvement']==(fu*up<oldL)
assert Q(result['tail_probability'])==tailprob
assert Q(result['weighted_box_width'])==fu*boxwidth
tail=fu*sum(Q(c['alpha'])for c in table['cones'])*496*pU**2*tailprob
assert Q(result['tail_upper'])==tail
assert Q(result['pi_factor_width'])==(fu-fl)*low
assert fu*up-fl*low==Q(result['weighted_box_width'])+Q(result['tail_upper'])+Q(result['pi_factor_width'])
assert Q(result['old128_upper_minus_new_upper'])==oldU-fu*up
assert Q(result['old_expression_lower_minus_new_upper'])==oldL-fu*up
assert fl*low<=oldU and fl*root_low<=oldU
assert result['global_lower_excludes_strict_saved_cap_pass']==(fl*low>=oldU)
assert result['early_stop_on_global_lower']==(fl*root_low>=oldU)
assert leafcounts==([1]*21 if result['early_stop_on_global_lower']else[256]*21)
assert result['old_beta2_interval']==spec['old_beta2_interval']
assert result['true_new_cap_is_only_upper']is True and result['family_failure_claim_from_nonpass']is False
assert not result['new_C6_bound']and not result['global_g0_proved']and not result['probability_renormalized']
print(json.dumps({'status':'PASS_SAVED_BOX_COVERAGE_POSITIVE_SUMS_AND_ARB768_KERNELS','leaves_checked':checked_boxes,'kernel_endpoint_checks':kernels,'precision_bits':768,'new_science_or_refinement':0}))
