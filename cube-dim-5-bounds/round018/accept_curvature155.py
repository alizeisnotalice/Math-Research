"""New complete g192 from audited153 margins and positive Taylor8; one C5 root."""
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
    a,c,old=read('audit157'),read('candidate155'),read('previous_upper')
    v=read('verified157')
    assert a['status']=='finite155_newq_log_jointcorners_and_complete110_cover_independently_qualified'
    assert a['objects']==107530 and a['leaves']==111201 and a['new_corners']==1779216
    assert a['old110_complete_domain_preserved'] and a['source_hashes_before_after_unchanged']
    assert c['complete_old110_domain_preserved'] and c['full110_leaf_count']==111201
    m=Q(a['minimum110'][0]);assert m==Q(c['minimum_margin110_new'])>0
    assert a['minimum110'][1]==c['minimum_witness_leaf_index']
    assert v['full110_closed_leaves']==111201 and Q(v['minimum_margin110_new'])==m
    assert v['leaf_index']==103387 and v['node']==[51742,'0111'] and v['epsilon']=='1/64'
    assert (a['passed115'],a['pending115'])==(109961,1240)
    assert (c['target115_passed'],c['target115_pending'])==(109961,1240)
    g0=Q(11,100000);term=Q(1);terms=[str(term)];poly=term
    for k in range(1,9):
        term=term*m/k;assert term>0;poly+=term;terms.append(str(term))
    exact_floor=g0*poly;N=1<<192
    units=(exact_floor.numerator*N)//exact_floor.denominator;g=Q(units,N)
    assert 0<=exact_floor-g<Q(1,N)
    linear=Q(c['candidate_floor']);assert linear==g0*(1+m)<g<=exact_floor
    I=Q(old['complete_I_upper']);oldg=Q(old['qualified_g']);old4=Q(old['best_upper']['fourth_power'])
    assert g>oldg>0 and I>0
    new4=3125*I**4/g**5;assert 0<new4<old4<32**4
    D=1<<256;k=isqrt(isqrt(new4.numerator*D**4//new4.denominator))
    lo,hi=Q(k,D),Q(k+1,D);pl,pu=map(Q,old['best_upper']['interval'])
    assert lo**4<=new4<=hi**4 and 0<lo<hi<pl
    simple=Q((hi.numerator*10**6+hi.denominator-1)//hi.denominator,10**6)
    assert hi<=simple<pl
    result={'status':'primary_accepts_complete153_g192_and_C5_upper','pid':os.getpid(),
        'seconds':time.monotonic()-t,'minimum_margin110':str(m),
        'minimum_witness_leaf_index':a['minimum110'][1],'minimum_witness_node':v['node'],
        'minimum_witness_epsilon':v['epsilon'],'positive_Taylor_degree':8,
        'Taylor_terms':terms,'g0':str(g0),'linear_floor':str(linear),
        'g_grid_bits':192,'g_grid_units':units,'qualified_g':str(g),
        'complete_I_upper':str(I),'exact_g_gain_vs_previous':str(g-oldg),
        'H128_floor_domain':'complete_F_only','G_floor_domain':'complete_Q_via_G_homogeneity',
        'best_upper':{'type':'unique_positive_fourth_root','fourth_power':str(new4),
                      'interval':[str(lo),str(hi)],'display_only':[float(lo),float(hi)]},
        'simple_rational_upper':str(simple),'exact_fourth_power_gain_vs_previous':str(old4-new4),
        'gain_interval_vs_previous':[str(pl-hi),str(pu-lo)],
        'comparisons':{n:{'old_upper':'32','exact_fourth_power_gain':str(Q(32)**4-new4),
                           'upper_gain_interval':[str(32-hi),str(32-lo)],'simple_exact_gain':str(32-simple)}
                       for n in ('checked_public','opening_internal','continuation_reference')},
        'closed_forest':{'roots':52360,'leaves':111201,'splits':58841,'pending':0},
        'target115':{'passed':109961,'pending':1240,'complete':False},
        'general_L1_all_positive_common_radii':True,'new_H_R_log_gradient_corner_or_exp_calls':0,
        'new_bound_root_isolations':1,'lower_changed':False,'dual_breakthrough_at':None,'post_dual_rounds':0,
        'proof':'round018/primary_curvature155_proof.md','sources':r['sources'],'registration_sha256':sha(reg)}
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'pid':result['pid'],'g_display_only':float(g),
                      'upper_display_only':result['best_upper']['display_only'],'simple_rational_upper':str(simple)}))
if __name__=='__main__':main()
