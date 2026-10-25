#!/usr/bin/env python3
"""Seven fixed pair-envelope checks, using independent binomial counts and asin."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from math import comb
import json,hashlib,os,time,sys
from flint import arb,fmpq,ctx
ctx.prec=640
sys.set_int_max_str_digits(0)
O=Path(__file__).resolve().parent
source=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round037/pair_resolved_green_envelope')
peaks=O.parents[1]/'round035/projection_peak_review/results.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,j):(O/n).write_text(json.dumps(j,indent=2)+'\n')
def B(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def saved_bounds(j):
 if 'common_binary_exponent'in j:
  e=F(2)**int(j['common_binary_exponent']);return F(int(j['lower_mantissa']))*e,F(int(j['upper_mantissa']))*e
 return tuple(F(int(j[w+'_mantissa']))*F(2)**int(j[w+'_binary_exponent']) for w in ['lower','upper'])
def overlap(x,j):
 a,b=map(F,ends(x));c,d=saved_bounds(j);assert max(a,c)<=min(b,d)
start=time.time();reg=json.loads((O/'input_manifest.json').read_text())
assert start<reg['shared_start_unix']+1200 and not(O/'execution_receipt.json').exists()
inputs=dict(reg['input_sha256'])
m=json.loads((source/'manifest.json').read_text())
for n,h in m['whitelist_sha256'].items():inputs[str(source/n)]=h
inputs.update(m['external_readonly_sha256'])
inputs[str(peaks)]=sha(peaks)
for p,h in inputs.items():assert sha(p)==h,p
assert sha(source/'check_diagnostic.py')=='e3cee920799d739b10be4ad54f3f2abacf1ef35f0b30346a7b5a5beb2d1dc78f'
r={'status':'RUNNING','pid':os.getpid(),'planned_configurations':7,'completed_configurations':0,'kernel_calls':0}
save('execution_receipt.json',r)
try:
 old=json.loads((source/'results.json').read_text());P={a['id']:F(a['R_exact']) for a in json.loads(peaks.read_text())['cases']}
 coef_expected=[{(4,1,1,0):10,(3,2,1,0):20},{(4,1,1,0):2,(3,2,1,0):8,(2,2,2,0):6,(3,1,1,1):8,(2,2,1,1):6},{(3,2,1,0):12,(2,2,1,1):18}]
 pi=arb.pi();M=[];rows=[]
 for k in (1,2,3):
  m=6-k;counts=Counter()
  for a in range(k+1):
   for b in range(m+1):
    if a in (0,k) and b in (0,m):continue
    typ=tuple(sorted((a,k-a,b,m-b),reverse=True))
    counts[typ]+=F(comb(k,a)*comb(m,b),2)
  assert dict(counts)==coef_expected[k-1] and sum(counts.values())==30
  s=sum((n*P['r3_'+''.join(map(str,p))] for p,n in counts.items()),F(0))
  aa=64*P[f'r2_k{k}'];cc=4096*s
  assert cc>960*aa
  # asin formula follows sin²(theta)=960A/C; coefficients are rational before pi^-6.
  value=(B(cc)/B(960).sqrt()*(B(960*aa/cc).sqrt()).asin()+B(aa*(cc-960*aa)).sqrt())/pi**6
  oo=old['configurations'][k-1]
  assert F(oo['S_exact'])==s and F(oo['A_times_pi6'])==aa and F(oo['C_times_pi6'])==cc
  overlap(value,oo['M']);M.append(value)
  rows.append({'configuration':k,'kind':'pair orbit','k':k,'third_class_counts':{str(p):str(n)for p,n in counts.items()},'S_exact':str(s),'A_times_pi6':str(aa),'C_times_pi6':str(cc),'M_bounds':ends(value),'display_only':str(value)})
  r['completed_configurations']+=1
 MT=192*M[0]+480*M[1]+320*M[2]
 assert 192+480+320==992;overlap(MT,old['global_M_ordered_sum'])
 for i,u in enumerate((F(1,4),F(1,2),F(3,4)),3):
  r2=3*(1+u*u);dv=F(3,2)*(1-u)**2;dw=F(3,2)*(1+u)**2
  coeff=1/r2*(1/dv+1/dw);assert coeff==F(4,9)/(1-u*u)**2
  lead=B(coeff)*M[2];collapse=MT/B(r2*dv);ratio=collapse/lead
  oo=old['configurations'][i]
  assert F(oo['u'])==u and F(oo['leading_coefficient_over_M3_exact'])==coeff
  overlap(lead,oo['pair_resolved_leading_coefficient']);overlap(collapse,oo['collapsed_common_distance_leading_coefficient']);overlap(ratio,oo['collapsed_over_resolved'])
  rows.append({'configuration':i+1,'kind':'fixed middle point','u':str(u),'exact_distance_squares':[str(r2),str(dv),str(dw)],'envelope_leading_over_M3_exact':str(coeff),'leading_bounds':ends(lead),'collapse_ratio_bounds':ends(ratio)})
  r['completed_configurations']+=1
 delta=F(1,64);q=F(6,5);v4=pi**2/2
 a=B(F(3,32)).sqrt()-B(3).sqrt()*B(delta);g=B(F(1,32)).sqrt()-B(3).sqrt()*B(delta)
 assert a>0 and g>0
 At=2*M[2]/a**2;BG=(MT-2*M[2])/(a**2*g)
 # Use radial integral (surface S3)/(4-q), then face cone/Jacobian/middle factor 1/2.
 radial=(4*v4)/(B(4-q))*B(delta)**B(4-q)/2
 volume=v4*B(delta)**4/2
 singular=At*radial**B(1/q);background=BG*volume**B(1/q);total=singular+background
 oo=old['configurations'][6]
 for value,key in [(a,'line_gap_a'),(g,'other_plane_gap_g'),(At,'tube_A_over_rho'),(BG,'background_B'),(singular,'singular_L6over5_upper'),(background,'background_L6over5_upper'),(total,'total_L6over5_upper')]:overlap(value,oo[key])
 assert total<B(F(1,11250)) and total<B(F(4,48750))
 rows.append({'configuration':7,'kind':'complete single middle tube','delta':'1/64','line_gap_bounds':ends(a),'other_plane_gap_bounds':ends(g),'singular_budget_bounds':ends(singular),'background_budget_bounds':ends(background),'total_budget_bounds':ends(total),'total_display_only':str(total),'over_old_target_bounds':ends(total*11250),'over_new_Bcrit_bounds':ends(total/B(F(4,48750))),'ordered_singular_terms':2,'ordered_background_terms':990})
 r['completed_configurations']+=1;assert r['completed_configurations']==len(rows)==7
 out={'status':'PASS_PAIR_RESOLVED_ENVELOPE_AND_COMPLETE_SINGLE_TUBE','configurations':rows,'q':'6/5','M_total_bounds':ends(MT),'whole_face_Lq_integrability':'proved analytically; norm not evaluated','true_G_asymptotic_claimed':False,'global_I_or_C6_certificate':False,'input_sha256':inputs,'source_sha256':sha(Path(__file__)),'precision_bits':ctx.prec,'pid':os.getpid(),'elapsed_seconds':time.time()-start,'kernel_calls':0,'old_enumerations':0,'grid':0,'Matlab':0}
 save('results.json',out);r['status']='PASS_STOPPED'
except BaseException as e:r['status']='FAILED_STOPPED_NO_RETRY';r['error']=repr(e);raise
finally:r['elapsed_seconds']=time.time()-start;save('execution_receipt.json',r)
print(json.dumps({'status':out['status'],'pid':out['pid'],'seconds':out['elapsed_seconds'],'tube':str(total),'singular':str(singular),'background':str(background),'new_target_ratio':str(total/B(F(4,48750)))},indent=2))
