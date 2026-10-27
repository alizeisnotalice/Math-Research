#!/usr/bin/env python3
"""N30B exact certificates and independent positive-quadrature diagnostic."""
import os,sys,time,json,hashlib,resource,concurrent.futures
from pathlib import Path
from fractions import Fraction as Q
sys.dont_write_bytecode=True
import mpmath as mp
from n30b import certify,ceilq
root=Path(__file__).resolve().parents[3]; search=json.loads(Path(__file__).with_name('n30b_search.json').read_text());table=json.loads((root/'revision_v3/numeric/table.json').read_text())
def run(row):
 start=time.monotonic();cpu=time.process_time();n=row['n']
 if row.get('status'):return row
 t=Q(row['t']);s=Q(row['sigma_candidate']);r=certify(n,t,s);old=Q(next(x['F_lower'] for x in table if x['n']==n));new=Q(r['new_F_upper'])
 U=None;precision=None
 for k in range(12,181,3):
  u=ceilq(new,10**k)
  if new<u<min(old,Q(r['old_F_lower'])):U=u;precision=k;break
 mp.mp.dps=80;mt=mp.mpf(t.numerator)/t.denominator;ms=mp.mpf(s.numerator)/s.denominator;R=mt*ms
 P=mp.quad(lambda x:n*(1-(1+1/(2*mt)-1/(2*x))**(n-1)/2)/x,[R,mt]);pl=Q(r['P_lower']);pu=Q(r['P_upper']);pm=mp.mpf(pl.numerator)/pl.denominator
 assert abs(P-pm)<mp.mpf('1e-55')
 r.update({'n':n,'t':str(t),'sigma':str(s),'R':row['R_candidate'],'original_table_F_lower':str(old),'strict_U':str(U) if U else None,'U_decimal_places':precision,'new_upper_lt_original_table_lower':new<old,'U_lt_original_table_lower':bool(U and U<old),'P_positive_integral_approx':mp.nstr(P,40),'positive_integral_check_abs_error_lt_1e-55':True,'gain_lower_decimal':mp.nstr(mp.mpf(Q(r['gain_lower']).numerator)/Q(r['gain_lower']).denominator,25),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
 if n==7:
  ref=json.loads((root/'upper/round028/cap_candidate/result.json').read_text()); m,e=ref['candidate_interval']['lower']; old28=Q(int(m))*Q(2)**int(e)
  assert U and U<old28
  r.update({'round028_full_lower_dyadic':[m,e],'strict_U_below_round028_full_lower':True})
 return r
if __name__=='__main__':
 start=time.monotonic()
 with concurrent.futures.ProcessPoolExecutor(max_workers=7) as ex: results=list(ex.map(run,search['results']))
 out={'batch':'N30B','stage':'exact certification','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'arithmetic_script_sha256':hashlib.sha256(Path(__file__).with_name('n30b.py').read_bytes()).hexdigest(),'search_script_sha256':search['script_sha256'],'wall_seconds':time.monotonic()-start,'workers':7,'results':results};Path(__file__).with_name('n30b_certificates.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps([{k:r.get(k) for k in ('n','P_minus_D_gt_tau_upper','new_upper_lt_old_lower','U_lt_original_table_lower','gain_lower_decimal','U_decimal_places','pid','wall_seconds')} for r in results],indent=2))
