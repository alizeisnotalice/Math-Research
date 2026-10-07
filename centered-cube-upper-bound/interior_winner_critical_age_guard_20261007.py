from pathlib import Path
from fractions import Fraction as F
import hashlib,json,time
P=Path('/Users/zhengzhihao/Desktop/T/output/cube_general_20261003/geom_20261006')
prefix='interior_winner_critical_age_guard_20261007'
reg={'purpose':'Full L1 quartic source, unique interior continuous winner; exact posterior age mean and half exponential moment, rational interval critical moment. No Lebesgue integration or nearmax claim.','source':'f=(1+epsilon/n sum(y_i^2-y_i^4))*1_[-2,2]^n','window':[1,2],'tau':'1/2','rounds':[[4,8,12],[16,64,256],[32,128,512]],'epsilon':['1/100','1/200'],'receivers':['all zero','all 1/10','alternating +/-1/10','first half 1/10, second half zero'],'arithmetic':'Fraction for mean and mgf at sigma1/2; log q interval 120 terms atanh tail; n multiple4','expected_tests':'mean<=1; half moment<=2; critical >= (352/401)*(1+n/2); q strictly between1 and4; all source density bounds kept','scope':'Formula verification for whole positive-Lebesgue interior family; selected discrete receivers check formulas only. Continuous theorem separate.','no_old_runs':True}
regpath=P/(prefix+'_registration.json')
assert not regpath.exists(),'registration already exists; refuse rerun'
regpath.write_text(json.dumps(reg,ensure_ascii=False,indent=2))
start=time.monotonic(); rows=[]
def log_interval(q):
    z=(q-1)/(q+1)
    assert 0<=z<1
    N=120
    lo=2*sum((z**(2*k+1)/F(2*k+1) for k in range(N)),F(0))
    return lo,lo+2*z**(2*N+1)/(F(2*N+1)*(1-z*z))
def sf(x):return str(x)
for ir, ns in enumerate(reg['rounds'],1):
  for n in ns:
    for eps in (F(1,100),F(1,200)):
      patterns=[[F(0)]*n,[F(1,10)]*n,[F((-1)**j,10) for j in range(n)],[F(1,10)]*(n//2)+[F(0)]*(n//2)]
      for ip,x in enumerate(patterns):
        s2=sum(t*t for t in x)/n;s4=sum(t**4 for t in x)/n
        C=1+eps*(s2-s4);A=eps*(F(1,12)-s2/2);B=eps/80
        q=A/(2*B);u=C+A*q-B*q*q
        def integral_power(k):return C*(q**(k//2)-1)/k+A*(q**((k+2)//2)-1)/(k+2)-B*(q**((k+4)//2)-1)/(k+4)
        mean=n*integral_power(n)/(q**(n//2)*u)
        half=1+F(n,2)*integral_power(n//2)/(q**(n//4)*u)
        loglo,loghi=log_interval(q)
        poly=A*(q-1)/2-B*(q*q-1)/4
        critlo=1+n*(C*loglo/2+poly)/u
        crithi=1+n*(C*loghi/2+poly)/u
        lo=1-12*eps;hi=1+eps/4
        lower=F(352,401)*(1+F(n,2))
        checks={'winner_internal':1<q<4,'q_range':F(47,15)<=q<=F(10,3),'source_density_lower':lo>=F(22,25),'source_density_upper':hi<=F(401,400),'weak_receiver':u>F(1,2),'mean_nonnegative':0<=mean,'mean_le_one':mean<=1,'half_between_one_two':1<=half<=2,'critical_linear_lower':critlo>=lower,'critical_le_one_B':crithi<=1+F(n,2)*loghi,'critical_interval_nonempty':critlo<=crithi}
        assert all(checks.values()),(ir,n,eps,ip,checks)
        rows.append({'round':ir,'n':n,'epsilon':sf(eps),'receiver_pattern':ip,'q':sf(q),'u_R':sf(u),'source_mass':'4^n * '+sf(1-F(28,15)*eps),'mean_exact':sf(mean),'half_mgf_exact':sf(half),'critical_lower_exact':sf(critlo),'critical_upper_exact':sf(crithi),'critical_linear_lower':sf(lower),'checks':checks,'diagnostic_floats':{'mean':float(mean),'half_mgf':float(half),'critical':float((critlo+crithi)/2),'critical_over_n':float((critlo+crithi)/(2*n)),'critical_interval_width':float(crithi-critlo)}})
result={'status':'PASS','rows':rows,'record_count':len(rows),'check_count':sum(len(r['checks']) for r in rows),'elapsed_seconds':time.monotonic()-start,'limitations':['No Monte Carlo','No source nearmax certification','No actual FIRST/history gate qualification','No weak-type constant lower bound; I/W may be tiny','Critical log rational enclosure is a certificate; floats are display only']}
respath=P/(prefix+'_results.json');respath.write_text(json.dumps(result,ensure_ascii=False,indent=2))
code=P/(prefix+'.py');code.write_bytes(Path(__file__).read_bytes())
receipt={'status':'terminal_exit_0','files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [regpath,respath,code]},'records':len(rows),'exact_checks':result['check_count'],'new_run':True}
(P/(prefix+'_receipt.json')).write_text(json.dumps(receipt,indent=2))
print(json.dumps({'status':'PASS','rows':len(rows),'checks':result['check_count'],'range_critical_over_n':[min(r['diagnostic_floats']['critical_over_n'] for r in rows),max(r['diagnostic_floats']['critical_over_n'] for r in rows)],'elapsed':result['elapsed_seconds']}))
