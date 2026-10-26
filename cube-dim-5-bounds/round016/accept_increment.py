"""One new primary acceptance and exact old/new bound comparison."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import hashlib,json,os,sys,time


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    regpath,out=map(Path,sys.argv[1:]);assert not out.exists()
    r=json.loads(regpath.read_bytes());start=time.monotonic()
    for name,row in r['sources'].items():assert sha(Path(row['path']))==row['sha256'],name
    def data(name):return json.loads(Path(r['sources'][name]['path']).read_bytes())
    old,new,audit=data('old_acceptance'),data('producer_result'),data('independent_result')
    g=Q(old['g']);I=Q(new['complete_I_upper']);oldI=Q(old['complete_I_upper'])
    assert I==Q(audit['complete_I_upper']) and 0<I<oldI
    cost4=3125*I**4/g**5;old4=Q(old['best_upper']['fourth_power'])
    assert cost4==Q(audit['new_cost_fourth_power_upper'])==Q(new['new_cost_fourth_power_upper'])<old4<32**4
    D=2**256;k=isqrt(isqrt(cost4.numerator*D**4//cost4.denominator))
    lo,hi=Q(k,D),Q(k+1,D);assert 0<lo<hi and lo**4<=cost4<=hi**4
    oldlo,oldhi=map(Q,old['best_upper']['interval']);assert hi<oldlo
    simple=Q((hi.numerator*10**6+hi.denominator-1)//hi.denominator,10**6)
    assert hi<=simple<oldlo
    result={'status':'primary_accepts_complete_C5_upper139','pid':os.getpid(),'seconds':time.monotonic()-start,
        'g':str(g),'complete_I_upper':str(I),'best_upper':{'type':'unique_positive_fourth_root','fourth_power':str(cost4),'interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)]},
        'simple_rational_upper':str(simple),'exact_I_gain_vs136':str(oldI-I),
        'exact_fourth_power_gain_vs136':str(old4-cost4),'upper_gain_interval_vs136':[str(oldlo-hi),str(oldhi-lo)],
        'comparisons':{name:{'old_upper':'32','new_upper_fourth_power':str(cost4),'strict_fourth_power_gain':str(Q(32)**4-cost4),'gain_interval':[str(32-hi),str(32-lo)],'simple_exact_gain':str(32-simple)} for name in ['checked_public','opening_internal','continuation_reference']},
        'same_general_L1_all_positive_common_radii':True,'g110_consumed':False,'lower_changed':False,'dual_breakthrough_at':None,'post_dual_rounds':0,
        'primary_proof':'round016/primary_increment_proof.md','root_bits_for_bound_isolation':256,'old_scalar_or_H_repeated':0,
        'registration_sha256':sha(regpath),'sources':r['sources']}
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'pid':result['pid'],'seconds':result['seconds'],'upper_display_only':result['best_upper']['display_only'],'simple_rational_upper':str(simple)}))


if __name__=='__main__':main()
