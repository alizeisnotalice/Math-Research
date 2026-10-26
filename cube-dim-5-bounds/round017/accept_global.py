"""One new complete g*/C5 acceptance; two exact bound isolations, no old H."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import hashlib,json,os,sys,time
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def root4(x):
    d=2**256;k=isqrt(isqrt(x.numerator*d**4//x.denominator));l,u=Q(k,d),Q(k+1,d)
    assert 0<l<u and l**4<=x<=u**4
    return l,u
def main():
    reg,out=map(Path,sys.argv[1:]);assert not out.exists();r=json.loads(reg.read_bytes());t=time.monotonic()
    for row in r['sources'].values():assert sha(Path(row['path']))==row['sha256']
    def read(n):return json.loads(Path(r['sources'][n]['path']).read_bytes())
    a,old=read('global148'),read('old139');assert a['pending']==0 and (a['passed'],a['splits'],a['roots'])==(111201,58841,52360)
    assert a['epsilon_leaf_counts']=={'1/128':12559,'1/64':98642}
    assert a['every_leaf_proves_same_epsilon_H_on_closed_box'] and a['every_admitted_rule_qualifies_saved_log_margin_as_whole_box_lower']
    assert a['all_source_hashes_unchanged'] and a['root_index_set_exact_and_all_roots_tiled']
    g0=Q(11,100000);m=Q(a['minimum_saved_log_margin']);g=g0*(1+m)
    assert m>0 and g==Q(a['optional_rational_G_cube_floor_for_primary'])==Q(a['optional_rational_face_H128_floor_for_primary'])
    assert len(a['minimum_margin_witnesses'])==1;w=a['minimum_margin_witnesses'][0];assert w['node']==[51752,'10']
    state=read('raw138');rows=[x for x in state['new_passed'] if x['node']==w['node']];assert len(rows)==1
    c=rows[0]['certificate'];assert digest(c)==w['certificate_sha256'] and c['passed'] and Q(c['minimum_log_margin'])==m
    assert c['target']=='11/100000' and c['leaf_epsilon']=='1/128'
    I=Q(old['complete_I_upper']);previous_g=Q(old['g']);assert g>g0>previous_g
    previous4=Q(old['best_upper']['fourth_power']);nominal4=3125*I**4/g0**5;best4=3125*I**4/g**5
    assert 0<best4<nominal4<previous4<32**4
    lo,hi=root4(best4);nl,nu=root4(nominal4);pl,pu=map(Q,old['best_upper']['interval'])
    assert hi<nl and nu<pl
    simple=Q((hi.numerator*10**6+hi.denominator-1)//hi.denominator,10**6);assert hi<=simple<nl
    result={'status':'primary_accepts_complete_global_gstar_and_C5_upper','pid':os.getpid(),'seconds':time.monotonic()-t,
        'nominal_g':str(g0),'qualified_g':str(g),'minimum_saved_log_margin':str(m),'minimum_witness':w,
        'H128_floor_domain':'complete_F_only','G_floor_domain':'complete_Q_via_G_homogeneity','complete_I_upper':str(I),
        'best_upper':{'type':'unique_positive_fourth_root','fourth_power':str(best4),'interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)]},
        'nominal_g110_cost':{'fourth_power':str(nominal4),'interval':[str(nl),str(nu)]},'simple_rational_upper':str(simple),
        'gain_interval_vs139':[str(pl-hi),str(pu-lo)],'exact_fourth_power_gain_vs139':str(previous4-best4),
        'gain_interval_vs_nominal_g110':[str(nl-hi),str(nu-lo)],
        'comparisons':{n:{'old_upper':'32','exact_fourth_power_gain':str(Q(32)**4-best4),'upper_gain_interval':[str(32-hi),str(32-lo)],'simple_exact_gain':str(32-simple)} for n in ['checked_public','opening_internal','continuation_reference']},
        'closed_forest':{'roots':52360,'passed':111201,'splits':58841,'pending':0},'epsilon_leaf_counts':a['epsilon_leaf_counts'],
        'general_L1_all_positive_common_radii':True,'new_H_PF_log_gradient_corner_or_exp_calls':0,'bound_root_isolations':2,
        'lower_changed':False,'dual_breakthrough_at':None,'post_dual_rounds':0,'proof':'round017/primary_global_proof.md',
        'source_sha256':r['sources'],'registration_sha256':sha(reg)}
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'pid':result['pid'],'seconds':result['seconds'],'g_display_only':float(g),'upper_display_only':result['best_upper']['display_only'],'simple_rational_upper':str(simple)}))
if __name__=='__main__':main()
