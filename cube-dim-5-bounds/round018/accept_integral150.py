"""One new complete I150/C5 isolation from frozen exact qualified inputs."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import hashlib, json, os, sys, time

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    reg, out = map(Path, sys.argv[1:])
    assert not out.exists()
    r = json.loads(reg.read_bytes()); t = time.monotonic()
    for name, row in r['sources'].items():
        assert sha(row['path']) == row['sha256'], name
    def read(name):
        return json.loads(Path(r['sources'][name]['path']).read_bytes())
    old, prod, audit = read('accepted017'), read('producer150'), read('audit152')
    final, inp, oldJ = read('final152'), read('inputs139'), read('result139')
    assert prod['pid'] == 69275 and audit['pid'] == 73479
    assert (prod['factors'],prod['roots'],prod['cached_band_products']) == (73,64,4672)
    assert (audit['factors'],audit['roots'],audit['cached_band_products']) == (73,64,4672)
    assert prod['all_sources_unchanged'] and audit['all_sources_before_after_unchanged']
    assert final['actual_audit_exitcode'] == 0 and final['actual_audit_pid_absent']
    before, after, gain = map(Q, (prod['old_UH'],prod['new_UH'],prod['exact_UH_gain']))
    assert before == Q(audit['old_UH']) == Q(inp['IH_upper'])
    assert after == Q(audit['new_UH']) and gain == Q(audit['exact_gain'])
    assert 0 < after < before and gain == before-after > 0
    oldI, g = Q(old['complete_I_upper']), Q(old['qualified_g'])
    J = Q(oldJ['auxiliary_J_interval'][1])
    assert oldI == Q(oldJ['complete_I_upper']) == before+J
    newI = after+J
    assert 0 < newI == oldI-gain < oldI
    old4 = Q(old['best_upper']['fourth_power'])
    new4 = 3125*newI**4/g**5
    assert 0 < new4 < old4 < Q(32)**4
    # A single new root; previous intervals are read, never recomputed.
    den = 1 << 256
    k = isqrt(isqrt(new4.numerator*den**4//new4.denominator))
    lo, hi = Q(k,den), Q(k+1,den)
    pl, pu = map(Q, old['best_upper']['interval'])
    assert 0 < lo < hi < pl and lo**4 <= new4 <= hi**4
    simple = Q((hi.numerator*10**6+hi.denominator-1)//hi.denominator, 10**6)
    assert hi <= simple < pl
    result = {
        'status':'primary_accepts_complete_I150_and_C5_upper',
        'pid':os.getpid(),'seconds':time.monotonic()-t,
        'qualified_g':str(g),'old_UH':str(before),'new_UH':str(after),
        'same_complete_J139_upper':str(J),'complete_I_upper':str(newI),
        'exact_I_gain_vs017':str(gain),
        'best_upper':{'type':'unique_positive_fourth_root','fourth_power':str(new4),
                      'interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)]},
        'simple_rational_upper':str(simple),
        'exact_fourth_power_gain_vs017':str(old4-new4),
        'gain_interval_vs017':[str(pl-hi),str(pu-lo)],
        'comparisons':{n:{'old_upper':'32','exact_fourth_power_gain':str(Q(32)**4-new4),
                          'upper_gain_interval':[str(32-hi),str(32-lo)],
                          'simple_exact_gain':str(32-simple)}
                       for n in ('checked_public','opening_internal','continuation_reference')},
        'general_L1_all_positive_common_radii':True,
        'H128_floor_domain':'complete_F_only','G_floor_domain':'complete_Q_via_G_homogeneity',
        'new_H_R_gradient_log_corner_or_old_power_calls':0,'new_bound_root_isolations':1,
        'lower_changed':False,'dual_breakthrough_at':None,'post_dual_rounds':0,
        'proof':'round018/primary_integral150_proof.md',
        'sources':r['sources'],'registration_sha256':sha(reg)}
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','pid','seconds','simple_rational_upper')} |
                     {'upper_display_only':result['best_upper']['display_only']}))

if __name__ == '__main__':
    main()
