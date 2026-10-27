"""One preregistered complete C5 bound; run only after 169 and 177 qualify."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import hashlib,json,os,sys,time

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    reg,out=map(Path,sys.argv[1:]); assert not out.exists()
    r=json.loads(reg.read_bytes()); start=time.monotonic()
    assert r['execution_ready'] is True
    assert r['first_audit_qualified'] is True
    assert r['code_sha256']==sha(__file__)
    for s in r['sources'].values(): assert sha(s['path'])==s['sha256']
    def read(k): return json.loads(Path(r['sources'][k]['path']).read_bytes())
    old,prod,audit=map(read,('accepted156','producer166','audit169'))
    # Exact schemas are bound in the final registration after the independent audit.
    assert audit['status']==r['audit169_success_status']
    assert audit[r['audit169_UH_field']]==prod['new_UH']
    assert audit['source_hashes_before_after_equal'] is True
    assert audit['worker_exitcodes']==[0]*15
    assert (audit['new_gain_cases'],audit['nodes'],audit['roots'])==(97971,97972,52360)
    assert audit['new_factor_mean_slope_witness_H_or_old_DFS_calls']==0
    assert r['audit169_result_qualification_reviewed'] is True
    assert prod['all_inputs_before_after_unchanged'] is True
    assert prod['new_gain_cases']==97971 and prod['nodes']==97972 and prod['roots']==52360
    assert prod['new_factor_mean_slope_H_R_log_oldpower_bandproduct_calls']==0
    g,J=Q(old['qualified_g']),Q(old['same_complete_J139_upper'])
    before,after=Q(old['new_UH']),Q(prod['new_UH'])
    assert 0<after<before and Q(prod['exact_gain_vs164'])>0
    newI=after+J; oldI=Q(old['complete_I_upper'])
    assert newI<oldI==before+J
    power=3125*newI**4/g**5; oldpower=Q(old['best_upper']['fourth_power'])
    assert 0<power<oldpower<32**4
    D=1<<256; k=isqrt(isqrt(power.numerator*D**4//power.denominator))
    lo,hi=Q(k,D),Q(k+1,D)
    assert lo**4<=power<hi**4
    oldlo,oldhi=map(Q,old['best_upper']['interval']); assert hi<oldlo
    simple=Q((hi.numerator*10**6+hi.denominator-1)//hi.denominator,10**6)
    assert hi<=simple<oldlo
    result={'status':'primary_accepts_complete_I166_and_C5_upper',
        'pid':os.getpid(),'seconds':time.monotonic()-start,
        'qualified_g':str(g),'new_UH':str(after),'same_complete_J139_upper':str(J),
        'complete_I_upper':str(newI),'exact_I_gain_vs156':str(before-after),
        'best_upper':{'type':'unique_positive_fourth_root','fourth_power':str(power),
            'interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)]},
        'simple_rational_upper':str(simple),'exact_fourth_power_gain_vs156':str(oldpower-power),
        'upper_gain_interval_vs156':[str(oldlo-hi),str(oldhi-lo)],
        'comparisons':{k:{'old_upper':'32','exact_fourth_power_gain':str(Q(32)**4-power),
            'upper_gain_interval':[str(32-hi),str(32-lo)],'simple_exact_gain':str(32-simple)}
            for k in ('checked_public','opening_internal','continuation_reference')},
        'general_L1_all_positive_common_radii':True,'H128_floor_domain':'complete_F_only',
        'G_floor_domain':'complete_Q_via_G_homogeneity','lower_changed':False,
        'dual_breakthrough_at':None,'post_dual_rounds':0,
        'new_H_R_gradient_log_corner_or_old_power_calls':0,'new_bound_root_isolations':1,
        'second_full_audit_status':'pending; this result is not the second audit',
        'sources':r['sources'],'registration_sha256':sha(reg)}
    with out.open('x') as f: json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'pid':result['pid'],
        'upper_display_only':result['best_upper']['display_only'],'simple':str(simple)}))

if __name__=='__main__': main()
