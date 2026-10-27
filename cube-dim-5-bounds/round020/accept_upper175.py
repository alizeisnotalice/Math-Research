"""One new full175 bound after independent175 audit and main proof review."""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import hashlib,json,os,sys,time

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 reg,out=map(Path,sys.argv[1:]);assert not out.exists();r=json.loads(reg.read_bytes());t=time.monotonic()
 assert r['execution_ready'] and r['main_theory_and_independent_result_reviewed']
 assert r['code_sha256']==sha(__file__)
 for source in r['sources'].values():assert sha(source['path'])==source['sha256']
 def read(k):return json.loads(Path(r['sources'][k]['path']).read_bytes())
 old,p,a=map(read,('accepted166','producer175','audit175'))
 assert p['status']=='candidate_pending_independent_audit_primary'
 assert a[r['audit_status_field']]==r['audit_success_status']
 assert a[r['audit_UH_field']]==p['new_UH']
 assert (a['new_fourth_roots'],a['new_slopes_gains'],a['nodes'],a['roots'])==(195942,97971,97972,52360)
 assert a['worker_exitcodes']==[0]*15 and a['all_source_hashes_before_after_equal'] is True
 assert p['all_inputs_before_after_unchanged'] and p['all_node_self_eff_monotone_vs166']
 assert p['all_new_root_left_nonnegative_next_positive']
 assert (p['new_fourth_roots'],p['new_gains'],p['nodes'],p['roots'])==(195942,97971,97972,52360)
 assert r['audit_real_new_roots_gains_full_tree_verified'] and r['all_actual_exitcodes_reviewed']
 UH=Q(p['new_UH']);oldUH=Q(old['new_UH']);assert UH<oldUH and Q(p['exact_gain_vs166'])==oldUH-UH
 g=Q(old['qualified_g']);J=Q(old['same_complete_J139_upper']);I=UH+J
 v=3125*I**4/g**5;oldv=Q(old['best_upper']['fourth_power']);assert v<oldv
 B=1<<256;n=isqrt(isqrt((v.numerator*B**4)//v.denominator));lo,hi=Q(n,B),Q(n+1,B)
 assert lo**4<=v<hi**4
 simple=Q((hi.numerator*10**6+hi.denominator-1)//hi.denominator,10**6);assert v<simple**4
 comparisons={}
 for label,old_upper in [('checked_public',Q(32)),('opening_internal',Q(32)),('continuation_U_ref',Q(32))]:
  assert v<old_upper**4;comparisons[label]={'old_upper':str(old_upper),'positive_fourth_power_gain':str(old_upper**4-v),'root_gain_interval':[str(old_upper-hi),str(old_upper-lo)]}
 result={'status':'primary_accepts_complete_I175_and_C5_upper','pid':os.getpid(),'seconds':time.monotonic()-t,'qualified_g':str(g),'new_UH':str(UH),'same_complete_J139_upper':str(J),'complete_I_upper':str(I),'exact_I_gain_vs166':str(oldUH-UH),'best_upper':{'type':'unique_positive_fourth_root','fourth_power':str(v),'interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)]},'simple_rational_strict_upper':str(simple),'exact_fourth_power_gain_vs166':str(oldv-v),'upper_gain_interval_vs166':[str(Q(old['best_upper']['interval'][0])-hi),str(Q(old['best_upper']['interval'][1])-lo)],'comparisons':comparisons,'general_L1_all_positive_common_radii':True,'lower_unchanged':True,'new_symbolic_or_numerical_old_H_R_J_g_log_calls':0,'source_registration_sha256':sha(reg)}
 out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','pid','seconds','simple_rational_strict_upper']}))
if __name__=='__main__':main()
