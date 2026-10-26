"""One new full I156/C5 certificate after finite161 suffix qualification."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import hashlib,json,os,sys,time
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    reg,out=map(Path,sys.argv[1:]);assert not out.exists()
    r=json.loads(reg.read_bytes());t=time.monotonic()
    for v in r['sources'].values():assert sha(v['path'])==v['sha256']
    def read(n):return json.loads(Path(r['sources'][n]['path']).read_bytes())
    old,prod,audit,a150=map(read,['previous155','producer156','audit161','accepted150'])
    assert prod['source_locks_before_after_unchanged'] and prod['pid']==77535
    assert audit['status']=='all292_factors_97907_means_and_fulltree_independently_qualified_pending_primary'
    assert audit['all_source_hashes_unchanged'] and audit['worker_exitcodes']==[0]*15
    assert audit['qualified_total_mean_checks']==97907 and audit['qualified_total_band_products']==7147211
    assert (audit['inherited_verified_factor_checks'],audit['inherited_verified_mean_checks'],audit['new_mean_checks'])==(292,41833,56074)
    assert audit['new_band_products']==4093402 and audit['new_factor_checks']==0
    assert audit['effective_improved_nodes']==prod['effective_improved_nodes']==97972
    before,after=Q(a150['new_UH']),Q(prod['new_UH'])
    assert Q(prod['accepted150_UH'])==before and Q(audit['new_UH'])==after
    gain=before-after;assert gain==Q(audit['exact_incremental_gain_vs150'])==Q(prod['exact_incremental_gain_vs150'])>0
    oldI,g=Q(old['complete_I_upper']),Q(old['qualified_g']);J=Q(a150['same_complete_J139_upper'])
    assert oldI==before+J
    newI=after+J;old4=Q(old['best_upper']['fourth_power']);new4=3125*newI**4/g**5
    assert 0<newI==oldI-gain<oldI and 0<new4<old4<32**4
    D=1<<256;k=isqrt(isqrt(new4.numerator*D**4//new4.denominator))
    lo,hi=Q(k,D),Q(k+1,D);pl,pu=map(Q,old['best_upper']['interval'])
    assert lo**4<=new4<=hi**4 and 0<lo<hi<pl
    simple=Q((hi.numerator*10**6+hi.denominator-1)//hi.denominator,10**6);assert hi<=simple<pl
    result={'status':'primary_accepts_complete_I156_and_C5_upper','pid':os.getpid(),'seconds':time.monotonic()-t,
        'qualified_g':str(g),'old_UH150':str(before),'new_UH':str(after),'same_complete_J139_upper':str(J),
        'complete_I_upper':str(newI),'exact_I_gain_vs155':str(gain),
        'best_upper':{'type':'unique_positive_fourth_root','fourth_power':str(new4),'interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)]},
        'simple_rational_upper':str(simple),'exact_fourth_power_gain_vs155':str(old4-new4),
        'gain_interval_vs155':[str(pl-hi),str(pu-lo)],
        'comparisons':{n:{'old_upper':'32','exact_fourth_power_gain':str(Q(32)**4-new4),
                          'upper_gain_interval':[str(32-hi),str(32-lo)],'simple_exact_gain':str(32-simple)}
                       for n in ('checked_public','opening_internal','continuation_reference')},
        'general_L1_all_positive_common_radii':True,'H128_floor_domain':'complete_F_only','G_floor_domain':'complete_Q_via_G_homogeneity',
        'new_H_R_gradient_log_corner_or_old_power_calls':0,'new_bound_root_isolations':1,
        'lower_changed':False,'dual_breakthrough_at':None,'post_dual_rounds':0,
        'proof':'round018/primary_integral156_proof.md','sources':r['sources'],'registration_sha256':sha(reg)}
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'pid':result['pid'],'upper_display_only':result['best_upper']['display_only'],'simple_rational_upper':str(simple)}))
if __name__=='__main__':main()
