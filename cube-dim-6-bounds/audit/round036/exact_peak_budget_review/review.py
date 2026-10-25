#!/usr/bin/env python3
"""Five fixed B036 crosschecks, zero kernel/enumeration/grid."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,os,time,sys
from flint import arb,fmpq,ctx
ctx.prec=640;sys.set_int_max_str_digits(0)
O=Path(__file__).resolve().parent
source=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round036/exact_peak_budget')
peaks=O.parents[1]/'round035/projection_peak_review/results.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,j):(O/n).write_text(json.dumps(j,indent=2)+'\n')
def B(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def oldends(row):return tuple(F(int(row[w+'_mantissa']))*F(2)**int(row[w+'_binary_exponent']) for w in ('lower','upper'))
def overlap(x,j):
 a,b=map(F,ends(x));c,d=oldends(j);assert max(a,c)<=min(b,d)
start=time.time();reg=json.loads((O/'input_manifest.json').read_text())
assert start<reg['shared_start_unix']+1200 and not(O/'execution_receipt.json').exists()
inputs=dict(reg['input_sha256']);m=json.loads((source/'manifest.json').read_text())
for n,h in m['whitelist_sha256'].items():inputs[str(source/n)]=h
inputs.update(m['external_readonly_sha256']);inputs[str(peaks)]=sha(peaks)
for p,h in inputs.items():assert sha(p)==h,p
r={'status':'RUNNING','pid':os.getpid(),'planned_configurations':5,'completed_configurations':0,'kernel_calls':0}
save('execution_receipt.json',r)
try:
 old=json.loads((source/'results.json').read_text());P={a['id']:F(a['R_exact'])for a in json.loads(peaks.read_text())['cases']}
 oo=old['configurations'];assert len(oo)==12
 pi=arb.pi();Cf=B(32*P['r1'])/pi**6;A=(Cf/pi**5)**B(F(1,6));K=B(F(6,7))*A
 overlap(Cf,oo[0]['Cf']);overlap(A,oo[0]['A']);overlap(K,oo[0]['K'])
 assert 7776*P['r1']==F(67320683,418119680)
 rows=[{'configuration':1,'kind':'Cf and interpolation coefficients','Cf_times_pi6':str(32*P['r1']),'Cf_bounds':ends(Cf),'A_bounds':ends(A),'K_bounds':ends(K),'Cf_new_over_old_exact':str(7776*P['r1'])}];r['completed_configurations']=1
 weights={'r3_4110':480,'r3_3210':1920,'r3_2220':480,'r3_3111':640,'r3_2211':1440}
 weighted=sum((n*P[k] for k,n in weights.items()),F(0));coef=24576*weighted
 assert weighted==F(44743079,66147840) and coef==F(715889264,43065)
 assert sum(weights.values())*6==29760
 B3=B(coef)/pi**6;overlap(B3,oo[6]['B3'])
 rows.append({'configuration':2,'kind':'complete ordered B3','weighted_R_exact':str(weighted),'B3_times_pi6':str(coef),'ordered_count':29760,'B3_bounds':ends(B3)});r['completed_configurations']=2
 eps=B(F(1,1000));tail=B(F(6,7))*(Cf*eps**7/pi**5)**B(F(1,6))
 overlap(tail,oo[8]['tail_upper_budget'])
 rows.append({'configuration':3,'kind':'whole-face short tail','epsilon':'1/1000','tail_bounds':ends(tail),'display_only':str(tail)});r['completed_configurations']=3
 goal=B(F(1,14000))+tail;cost=(6*goal*65000)**B(F(6,5));Igoal=goal**B(F(6,5));crit=F(4,48750)
 assert goal<B(crit) and cost<64
 assert old['new_global_g_lower_bound']is False and old['new_global_H_norm_upper_bound']is False
 for x,k in [(goal,'B_goal'),(cost,'conditional_cost'),(Igoal,'I_goal')]:overlap(x,oo[10][k])
 rows.append({'configuration':4,'kind':'conditional cost only','g_target':'1/65000','H_norm_goal':'1/14000','B_critical_exact':str(crit),'B_goal_bounds':ends(goal),'I_goal_bounds':ends(Igoal),'conditional_cost_bounds':ends(cost),'cost_display_only':str(cost),'both_global_premises_proved':False});r['completed_configurations']=4
 delta=B(F(1,16));V=pi**2/2*(10+10*B(2).sqrt()+2*B(5).sqrt())*delta**4
 overlap(V,oo[7]['volume_upper']);e=(Cf*pi*V)**B(F(1,5))
 # Positive arctangent complementary-angle identity, independent of pi/2-atan(32e).
 R=1/(3*e**3)-1024/e+32768*(1/(32*e)).atan()
 short=K*e**B(F(7,6));long=V**B(F(5,6))*Cf*R
 U=short+long;L=short/(1+1024*e**2)+long
 for x,k in [(e,'epsilon_star'),(R,'R'),(U,'allocated_upper_budget'),(L,'best_two_envelope_lower_budget')]:overlap(x,oo[11][k])
 assert L/B(crit)>B(F(237041,62500))>1
 rows.append({'configuration':5,'kind':'unique fixed tube budget diagnostic','delta':'1/16','volume_bounds':ends(V),'epsilon_rule_bounds':ends(e),'allocated_upper_bounds':ends(U),'best_two_envelope_lower_bounds':ends(L),'method_lower_ratio_bounds':ends(L/B(crit)),'strict_safe_method_ratio':'237041/62500','actual_G_norm_lower_bound':False});r['completed_configurations']=5
 assert len(rows)==5
 out={'status':'PASS_EXACT_CONSTANTS_TAIL_CONDITIONAL_COST_AND_CAP_DIAGNOSTIC','configurations':rows,'input_sha256':inputs,'source_sha256':sha(Path(__file__)),'precision_bits':ctx.prec,'pid':os.getpid(),'elapsed_seconds':time.time()-start,'new_global_g':False,'new_global_H_norm':False,'new_I_or_C6_certificate':False,'off_S2_numeric_configuration_replayed':False,'kernel_calls':0,'old_enumerations':0,'grid':0,'Matlab':0}
 save('results.json',out);r['status']='PASS_STOPPED'
except BaseException as e:r['status']='FAILED_STOPPED_NO_RETRY';r['error']=repr(e);raise
finally:r['elapsed_seconds']=time.time()-start;save('execution_receipt.json',r)
print(json.dumps({'status':out['status'],'pid':out['pid'],'seconds':out['elapsed_seconds'],'Cf':str(Cf),'B3':str(B3),'tail':str(tail),'conditional_cost':str(cost),'cap_lower_ratio':str(L/B(crit))},indent=2))
