#!/usr/bin/env python3
"""O31A: one exact candidate, full R030 interval comparison and finite-a check."""
import os,sys,json,time,hashlib,resource
from fractions import Fraction as Q
from pathlib import Path
sys.dont_write_bytecode=True
path=Path(__file__).resolve();clone=path.parents[3];sys.path.insert(0,str(clone/'upper/round030/compute'));sys.path.insert(0,str(clone/'revision_v3/numeric'))
from n30b import certify,Pbounds,Bbounds,sqrtbounds,ceilq
from verify import log_bounds
started=time.monotonic();cpu=time.process_time();search=json.loads(path.with_name('o31a_search.json').read_text());found=[r for r in search['results'] if r['best']];row=min(found,key=lambda x:float(x['best']['F_float']));b=row['best'];n=7;t=Q(b['t_rational']);s=Q(b['sigma_rational']);R=t*s
oldpath=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/audit/round030/numeric/independent_n7.json');gold=json.loads(oldpath.read_text());oldlo=Q(gold['F_new'][0]);result=certify(n,t,s)
assert result['P_minus_D_gt_tau_upper'] and result['new_upper_lt_old_lower'];newhi=Q(result['new_F_upper']);U=None
for k in range(12,181,3):
 q=ceilq(newhi,10**k)
 if newhi<q<oldlo:U=q;break
assert U is not None
# Scalar finite-a bound at one rational guard.
a=2**512;u=Q(1,a);alpha=Q(a,a+1);beta=1+1/t;Rg=R+Q(1,2**80)
assert Rg<t and alpha*beta>1+s/Rg+Q(2,a+1)
pl,_=Pbounds(n,t,Rg/t);Bl,Bu=Bbounds(n,t);D=Q(result['D']);logu=log_bounds(t)[1]
Eup=Bu-D+((1+u)**n-1)*(1+n*logu)+u**n;Aa=(beta+u)**n;tau=sqrtbounds(Aa*Eup)[1];margin=alpha**(-n)*pl-tau-alpha**(-n)*D
assert margin>0
finite=(1+Aa/tau)*(tau+Eup);rounded=ceilq(finite,10**150);assert rounded<U<oldlo
out={'batch':'O31A','chosen_cell':row['cell'],'t':str(t),'R_t_sigma':str(R),'sigma':str(s),'g_positive_by_exact_intervals':True,'B_lower':result['B_lower'],'B_upper':result['B_upper'],'D':result['D'],'P_lower':result['P_lower'],'P_upper':result['P_upper'],'tau_upper':result['tau_upper'],'F_lower':result['new_F_lower'],'F_upper':result['new_F_upper'],'U':str(U),'round030_F_new_full_lower':str(oldlo),'strict_U_below_round030_full_lower':True,'strict_gap_R030_lower_minus_new_upper':str(oldlo-newhi),'gap_decimal':format(float(oldlo-newhi),'.16g'),'finite_a':'2^512','finite_guard_R':str(Rg),'finite_guard_positive':True,'finite_budget_margin_positive':True,'finite_budget_margin_decimal':format(float(margin),'.15g'),'finite_upper_rounded':str(rounded),'finite_upper_below_U':True,'global_optimality_claim':False,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'search_sha256':hashlib.sha256(path.with_name('o31a_search.py').read_bytes()).hexdigest(),'arithmetic_sha256':hashlib.sha256((clone/'upper/round030/compute/n30b.py').read_bytes()).hexdigest(),'old_independent_certificate_sha256':hashlib.sha256(oldpath.read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-started,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
path.with_name('o31a_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k in ['chosen_cell','t','R_t_sigma','sigma','g_positive_by_exact_intervals','U','strict_U_below_round030_full_lower','gap_decimal','finite_guard_positive','finite_budget_margin_positive','finite_budget_margin_decimal','finite_upper_below_U','pid','wall_seconds']},indent=2))
