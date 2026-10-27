#!/usr/bin/env python3
"""One fixed finite-a guard and exact scalar certificate for each N30B point."""
import sys,os,time,json,hashlib,resource,concurrent.futures
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True
from n30b import Pbounds,Bbounds,sqrtbounds,ceilq
from verify import log_bounds
path=Path(__file__).resolve();inp=json.loads(path.with_name('n30b_certificates.json').read_text());ROWS={r['n']:r for r in inp['results']}
def run(n):
 start=time.monotonic();cpu=time.process_time();row=ROWS[n];t=Q(row['t']);s=Q(row['sigma']);U=Q(row['strict_U']);a=2**512;u=Q(1,a);alpha=Q(a,a+1);beta=1+1/t;D=Q(row['D']);R=t*s+Q(1,2**80)
 assert t*s<R<t
 guard=alpha*beta-(1+s/R+Q(2,a+1));assert guard>0
 pl,_=Pbounds(n,t,R/t);Bl,Bu=Bbounds(n,t);logu=log_bounds(t)[1]
 delta=((1+u)**n-1)*(1+n*logu)+u**n
 Eup=Bu-D+delta;Aa=(beta+u)**n;assert Eup>0
 tau=sqrtbounds(Aa*Eup)[1]; plfin=alpha**(-n)*pl; df=alpha**(-n)*D
 condition=plfin-tau-df; assert condition>0
 finite=(1+Aa/tau)*(tau+Eup);assert finite<U
 # Round outward to limit saved data size; preserve strict comparisons.
 rounded=ceilq(finite,10**160);assert rounded<U
 return {'n':n,'pid':os.getpid(),'t':str(t),'sigma':str(s),'a':'2^512','R_guard':str(R),'R_minus_t_sigma':'1/2^80','strict_guard_positive':True,'finite_P_lower_minus_tau_minus_scaled_D_positive':True,'finite_cost_upper_rounded':str(rounded),'strict_U':str(U),'finite_cost_upper_lt_U':True,'strict_U_lt_original_table_F_lower':row['U_lt_original_table_lower'],'guard_margin_decimal':format(float(guard),'.12g'),'budget_margin_decimal':format(float(condition),'.12g'),'U_minus_finite_decimal':format(float(U-finite),'.12g'),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
if __name__=='__main__':
 dims=[7] if len(sys.argv)==1 else [4,5,6,8,16,64];start=time.monotonic()
 if len(dims)==1:res=[run(7)]
 else:
  with concurrent.futures.ProcessPoolExecutor(max_workers=6) as ex:res=list(ex.map(run,dims))
 out={'batch':'N30B finite','workers':len(dims),'script_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'input_certificates_sha256':hashlib.sha256(path.with_name('n30b_certificates.json').read_bytes()).hexdigest(),'wall_seconds':time.monotonic()-start,'results':res}
 suffix='n7' if len(dims)==1 else 'other';path.with_name('n30b_finite_'+suffix+'.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(res,indent=2))
