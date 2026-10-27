#!/usr/bin/env python3
"""Primary integration: recheck the seven frozen certificates, no optimization."""
import sys,json,time,hashlib,os,resource,concurrent.futures
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent/'compute'))
from n30b import certify,sqrtbounds,ceilq
INPUT=Path(__file__).resolve().parent/'compute/n30b_certificates.json'
ROWS=json.loads(INPUT.read_text())['results']
def check(row):
 start=time.monotonic();cpu=time.process_time();n=row['n'];t=Q(row['t']);s=Q(row['sigma']);U=Q(row['strict_U'])
 x=certify(n,t,s)
 for k in ['B_lower','B_upper','P_lower','P_upper','D','tau_upper','new_F_lower','new_F_upper','old_F_lower','old_F_upper']:
  assert x[k]==row[k],(n,k)
 assert x['P_minus_D_gt_tau_upper'] and x['new_upper_lt_old_lower']
 assert Q(x['new_F_upper'])<U<Q(row['original_table_F_lower'])
 # Independently bound the finite-a tail loss by its positive integrand.
 a=2**512; u=Q(1,a);alpha=Q(a,a+1);beta=1+1/t;D=Q(x['D']);R0=t*s;R=R0+Q(1,2**80)
 assert 1<R0<R<t
 guard=alpha*beta-1-s/R-Q(2,a+1);assert guard>0
 Pguard=Q(x['P_lower'])-n*(R-R0)/R0
 error=((1+u)**n-1)*(1+n*(t-1))+u**n
 E=Q(x['B_upper'])-D+error;Aa=(beta+u)**n;tau=sqrtbounds(Aa*E)[1]
 assert alpha**(-n)*(Pguard-D)>tau
 finite=(1+Aa/tau)*(tau+E);assert finite<U
 if n==7:
  raw=json.loads((ROOT/'upper/round028/cap_candidate/result.json').read_text())['candidate_interval']['lower']
  assert U<Q(int(raw[0]))*Q(2)**int(raw[1])
 return {'n':n,'strict_U':str(U),'fixed_point_and_finite_a':'PASS','finite_cost_upper':str(ceilq(finite,10**80)),'a':'2^512','R_increment':'2^-80','pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
if __name__=='__main__':
 start=time.monotonic()
 with concurrent.futures.ProcessPoolExecutor(max_workers=7) as e: results=list(e.map(check,ROWS))
 sources=[Path(__file__),INPUT,INPUT.with_name('n30b.py')]
 out={'status':'PASS primary integration, seven accepted points and finite-a certificates','no_new_search':True,'workers':7,'wall_seconds':time.monotonic()-start,'results':results,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
 Path(__file__).with_name('primary_verification.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASS seven accepted points, finite a=2^512, complete old-cost comparisons; no search')
