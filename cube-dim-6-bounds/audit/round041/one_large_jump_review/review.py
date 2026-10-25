#!/usr/bin/env python3
"""F041: five independent fixed checks, no old frame enumeration or spatial kernel."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json,hashlib,os,time,sys
import sympy as sp
from flint import arb,fmpq,ctx
ctx.prec=640;sys.set_int_max_str_digits(0)
O=Path(__file__).resolve().parent
source=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round041/one_large_jump_cap')
peaks=O.parents[1]/'round035/projection_peak_review/results.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,j):(O/n).write_text(json.dumps(j,indent=2)+'\n')
def B(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def oldends(j):return tuple(F(int(j[w+'_mantissa']))*F(2)**int(j[w+'_binary_exponent']) for w in ('lower','upper'))
def overlap(x,j):
 a,b=map(F,ends(x));c,d=oldends(j);assert max(a,c)<=min(b,d)
start=time.time();reg=json.loads((O/'input_manifest.json').read_text())
assert start<reg['started_unix']+900 and not(O/'execution_receipt.json').exists()
inputs=dict(reg['input_sha256']);m=json.loads((source/'manifest.json').read_text())
for n,h in m['whitelist_sha256'].items():inputs[str(source/n)]=h
inputs.update(m['external_readonly_sha256'])
inputs[str(peaks)]='9a960d07a5fba0932788008a37e77516536d03a62c33eb1ddc946b1cd5f0f647'
for p,h in inputs.items():assert sha(p)==h,p
r={'status':'RUNNING','pid':os.getpid(),'planned_configurations':5,'completed_configurations':0,'old_frame_enumerations':0,'G_H_native_calls':0}
save('execution_receipt.json',r)
try:
 old=json.loads((source/'results.json').read_text());oo=old['configurations'];assert len(oo)==6
 # 1. Symbolic all-pair formula, with one exact minimizing witness.
 n0,n1,n2,n3=sp.symbols('n0 n1 n2 n3');N=n0+n1+n2+n3
 c01=n0+n1-n2-n3;c02=n0-n1+n2-n3;c12=n0-n1-n2+n3
 Gram=sp.Matrix([[N,c01,c02],[c01,N,c12],[c02,c12,N]])
 e3=n0*n1*n2+n0*n1*n3+n0*n2*n3+n1*n2*n3
 aa=n0+n3;bb=n1+n2
 assert sp.expand(Gram.det()-16*e3)==0
 assert sp.expand(N*N-c12*c12-4*aa*bb)==0
 assert sp.cancel(Gram.det()/(N*N-c12*c12)-4*(n0*n3/aa+n1*n2/bb))==0
 h=sp.symbols('h');assert sp.expand(h*(6-h)-5-(h-1)*(5-h))==0
 witness=sp.Matrix([[6,4,4],[4,6,2],[4,2,6]])
 assert witness.det()==64 and sp.Matrix([[6,2],[2,6]]).det()==32
 assert sp.Rational(6)-sp.Rational(16,6)==sp.Rational(10,3)
 rows=[{'configuration':1,'kind':'global separation proof and one equality witness','line_min_squared':'10/3','pair_min_squared':'2','all_pair_Gram_identity':'distance^2=4*(n0*n3/(n0+n3)+n1*n2/(n1+n2))','witness_triple_Gram':64,'witness_pair_Gram':32,'coverage':'all31 lines/all465 pairs analytically; no frame enumeration'}];r['completed_configurations']=1
 assert oo[0]['line_count']==31 and oo[0]['pair_count']==465
 assert oo[0]['line_min_squared']=='10/3' and oo[0]['pair_min_squared']=='2'
 # 2. Transitivity incidence counts, no 4495-triple replay.
 data=json.loads(peaks.read_text());P={x['id']:x for x in data['cases']}
 names=['4110','3210','2220','3111','2211'];restricted={};weighted=F(0);allweighted=F(0)
 for name in names:
  case=P['r3_'+name];nn=case['triple_orbit_count'];rr=F(case['R_exact'])
  np=F(29,32)*nn;assert np.denominator==1
  restricted[name]=np.numerator;weighted+=np*rr;allweighted+=nn*rr
 assert restricted==dict(zip(names,[435,1740,435,580,1305]))
 assert sum(restricted.values())==comb(31,3)==4495 and 6*sum(restricted.values())==31*30*29==26970
 coef=24576*weighted;assert coef==F(44743079,2970)
 assert coef==F(29,32)*F(715889264,43065)
 assert restricted==oo[1]['unordered_counts'] and coef==F(oo[1]['Crem_times_pi6_exact'])
 rows.append({'configuration':2,'kind':'exact restricted triple incidence','unordered_counts':restricted,'unordered_total':4495,'ordered_total':26970,'Crem_times_pi6':str(coef),'ratio_to_B3':'29/32'});r['completed_configurations']=2
 # 3. One selected-parameter coefficient, no new quotient peak computation.
 pi=arb.pi();R1=F(P['r1']['R_exact']);D=B(R1/32)/pi**6;Cbig=B(R1/8)/pi**6
 assert (R1/8)/(32*R1)==F(1,256)
 overlap(Cbig,oo[3]['Cbig'])
 rows.append({'configuration':3,'kind':'one large parameter coefficient','D_times_pi6':str(R1/32),'Cbig_times_pi6':str(R1/8),'Cbig_over_Cf_exact':'1/256','D_bounds':ends(D),'Cbig_bounds':ends(Cbig)});r['completed_configurations']=3
 # 4. The sole rho cap; alternative factored c0 and squared-volume expression.
 rho=F(1,64);eta=B(6).sqrt()*B(rho)
 c0=B(6).sqrt()*B(F(1,2)-rho);c1=B(F(5,6)).sqrt()-eta;c2=B(F(1,2)).sqrt()-eta
 vol=(B(6).sqrt()/6)*(8*pi**2/15)*B(rho)**5
 assert c0>0 and c1>0 and c2>0
 for value,key in [(eta,'eta'),(c0,'c0'),(c1,'c1'),(c2,'c2'),(vol,'cap_volume')]:overlap(value,oo[2][key])
 assert oo[2]['time_denominator_coefficients']==[992,960,928]
 rows.append({'configuration':4,'kind':'sole cap geometry constants','rho':'1/64','eta_bounds':ends(eta),'c0_bounds':ends(c0),'c1_bounds':ends(c1),'c2_bounds':ends(c2),'cap_volume_bounds':ends(vol),'time_denominator_coefficients':[992,960,928],'frame_lower_eigenvalue':26,'frame_upper_eigenvalue':32});r['completed_configurations']=4
 # 5. One original shortcap epsilon, all complementary events paid.
 eps=F(1,1000);Crem=B(coef)/pi**6;Mrem=Crem/(c0*c1*c2)**2
 good=B(F(6,7))*(Cbig*B(eps)**7/pi**5)**B(F(1,6))
 bad=B(eps)*Mrem*(vol**5)**B(F(1,6));total=good+bad
 overlap(good,oo[3]['good_branch_short_norm_upper']);overlap(Crem,oo[4]['Crem']);overlap(Mrem,oo[4]['uniform_bad_density_upper']);overlap(bad,oo[4]['bad_branch_short_norm_upper']);overlap(total,oo[5]['short_cap_norm_upper'])
 simple=F(1948138674853,500000000000000000)
 assert F(ends(total)[1])<simple==F(oo[5]['simple_rational_upper'])
 assert total<B(F(4,48750))
 assert oo[5]['time_tail_above_epsilon_paid']is False and oo[5]['entire_face_integral_claimed']is False
 rows.append({'configuration':5,'kind':'single original short potential cap','epsilon':'1/1000','threshold_abs_a':'1/2','Crem_bounds':ends(Crem),'Mrem_bounds':ends(Mrem),'good_branch_bounds':ends(good),'entire_complement_bounds':ends(bad),'total_shortcap_bounds':ends(total),'strict_simple_rational_upper':str(simple),'over_conditional_Bcrit_bounds':ends(total/B(F(4,48750))),'time_tail_paid':False});r['completed_configurations']=5
 assert len(rows)==5
 result={'status':'PASS_ALL_TIME_POINTWISE_TWO_BRANCH_BOUND_AND_ONE_SHORTCAP','configurations':rows,'input_sha256':inputs,'source_sha256':sha(Path(__file__)),'precision_bits':ctx.prec,'pid':os.getpid(),'elapsed_seconds':time.time()-start,'all_positive_t_pointwise_inequality':'proved analytically including pointwise_coarea_note','shortcap_only_epsilon':'1/1000','full_G_cap_certified':False,'global_I_g_C6_certified':False,'new_G_H_native_grid_old_enumeration_MATLAB':0}
 save('results.json',result);r['status']='PASS_STOPPED'
except BaseException as e:r['status']='FAILED_STOPPED_NO_RETRY';r['error']=repr(e);raise
finally:r['elapsed_seconds']=time.time()-start;save('execution_receipt.json',r)
print(json.dumps({'status':result['status'],'pid':result['pid'],'seconds':result['elapsed_seconds'],'good':str(good),'bad':str(bad),'total':str(total),'Mrem':str(Mrem),'simple_upper':str(simple)},indent=2))
