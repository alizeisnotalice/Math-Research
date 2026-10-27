"""One final175 audit: integer cross-products only; no root/H/node replay."""
from pathlib import Path
from math import gcd
from datetime import datetime,timezone
import hashlib,json,os,sys,time
def read(p):return json.loads(Path(p).read_bytes())
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()
def q(s):
    z=s.split('/');assert len(z) in (1,2)
    a,b=int(z[0]),int(z[1]) if len(z)==2 else 1
    assert b>0;return a,b
def equal(a,b):return a[0]*b[1]==b[0]*a[1]
def less(a,b):return a[0]*b[1]<b[0]*a[1]
def add(a,b):return a[0]*b[1]+b[0]*a[1],a[1]*b[1]
def sub(a,b):return a[0]*b[1]-b[0]*a[1],a[1]*b[1]
def textq(a):
    x,y=a;g=gcd(x,y);x//=g;y//=g
    return str(x) if y==1 else f'{x}/{y}'
def main():
    regpath,out=map(Path,sys.argv[1:]);assert not out.exists()
    reg=read(regpath);assert reg['maximum_math_runs']==1
    assert reg['method']=='integer_cross_products_no_Fraction_SymPy_or_root_algorithm'
    assert sha(__file__)==reg['code_sha256']
    manifest=read(reg['input_manifest']['path']);assert sha(reg['input_manifest']['path'])==reg['input_manifest']['sha256']
    start=time.monotonic();before={p:sha(p) for p in manifest['files']}
    assert all(before[p]==v['sha256'] for p,v in manifest['files'].items())
    d,old,prod,audit,mat,mi=(read(reg['named_inputs'][k]) for k in ['accepted175','accepted166','producer175','audit175','matlab_result','matlab_inputs'])
    checks={}
    def check(k,value):assert value,k;checks[k]=True
    check('formal_exact175_status',d['status']=='primary_accepts_complete_I175_and_C5_upper')
    check('actual_component_identity',d['new_UH']==prod['new_UH']==audit['new_UH'])
    check('same_g192_no173_consumed',d['qualified_g']==old['qualified_g'])
    check('same_complete_J139',d['same_complete_J139_upper']==old['same_complete_J139_upper'])
    check('audit_counts_and_exitcodes',(audit['new_fourth_roots'],audit['new_slopes_gains'],audit['nodes'],audit['roots'])==(195942,97971,97972,52360) and audit['worker_exitcodes']==[0]*15)
    check('audited_source_identity',audit['all_source_hashes_before_after_equal'] and prod['all_inputs_before_after_unchanged'])
    U,J,I,g=(q(d[k]) for k in ['new_UH','same_complete_J139_upper','complete_I_upper','qualified_g'])
    check('I_equals_full_UH_plus_J',equal(I,add(U,J)))
    check('positive_complete_I_and_g',I[0]>0 and g[0]>0 and U[0]>0 and J[0]>0)
    check('new_I_strictly_below166',less(I,q(old['complete_I_upper'])))
    check('exact_I_gain_matches_UH_difference',equal(q(d['exact_I_gain_vs166']),sub(q(old['new_UH']),U)))
    # Clear denominators once.  This differs from the accepted Fraction/isqrt path.
    V=q(d['best_upper']['fourth_power']);In,Id=I;gn,gd=g
    expected=(3125*In**4*gd**5,Id**4*gn**5)
    check('complete_fourth_power_identity',equal(V,expected))
    check('strict_fourth_power_gain_vs166',less(V,q(old['best_upper']['fourth_power'])))
    check('saved_exact_fourth_power_gain',equal(q(d['exact_fourth_power_gain_vs166']),sub(q(old['best_upper']['fourth_power']),V)))
    lo,hi=map(q,d['best_upper']['interval']);check('positive_ordered_saved_interval',0<lo[0] and less(lo,hi))
    vl=V[0]*lo[1]**4-lo[0]**4*V[1]
    vu=hi[0]**4*V[1]-V[0]*hi[1]**4
    check('saved_root_lower_by_cross_power',vl>=0)
    check('saved_root_upper_by_cross_power',vu>0)
    oldlo,oldhi=map(q,old['best_upper']['interval'])
    check('new_saved_interval_strictly_below_old',less(hi,oldlo))
    for k,want in zip(d['upper_gain_interval_vs166'],[sub(oldlo,hi),sub(oldhi,lo)]):
        check('root_gain_endpoint_'+str(len(checks)),equal(q(k),want) and want[0]>0)
    C=q(d['simple_rational_strict_upper'])
    check('stated_portable_is_exact',equal(C,(14460339,500000)))
    check('portable_strict_by_cross_power',C[0]**4*V[1]-V[0]*C[1]**4>0)
    check('portable_covers_saved_upper',not less(C,hi))
    for name,row in d['comparisons'].items():
        O=q(row['old_upper']);gain=(O[0]**4*V[1]-V[0]*O[1]**4,O[1]**4*V[1])
        check(name+'_correct_old32',equal(O,(32,1)))
        check(name+'_positive_exact4gain',gain[0]>0 and equal(gain,q(row['positive_fourth_power_gain'])))
        for i,want in enumerate([sub(O,hi),sub(O,lo)]):check(name+'_rootgain'+str(i),equal(q(row['root_gain_interval'][i]),want))
    check('MATLAB_actual_eight_true',mat['status']=='PASS_8_new175_exact_scalar_checks' and mat['checks']==[True]*8)
    check('MATLAB_same_candidate',mat['candidate_input_sha256']==sha(reg['named_inputs']['producer175'])==mi['candidate_sha256'])
    check('MATLAB_exact_UH_and_gain_binding',equal(q(mi['candidate_UH']),U) and equal(q(mi['gain']),q(d['exact_I_gain_vs166'])))
    check('formal_scope_preserved',d['general_L1_all_positive_common_radii'] and d['lower_unchanged'])
    after={p:sha(p) for p in manifest['files']};check('all_frozen_sources_unchanged',before==after)
    result={'status':'PASS_frozen175_final_exact_upper_audit','pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),
      'seconds':time.monotonic()-start,'checks':checks,'input_files':len(before),'all_source_locks_equal_before_after':True,
      'complete_I':d['complete_I_upper'],'qualified_g':d['qualified_g'],'complete_UH':d['new_UH'],
      'portable_strict_C5_upper':d['simple_rational_strict_upper'],
      'saved_interval_cross_power_residuals':{'lower':str(vl),'upper':str(vu)},
      'root_search_or_isolation_calls':0,'old_H_R_J_mean_g_slope_nodegain_DFS_baseline_calls':0,
      'new_current_final_constant_crosscheck_runs':1,'173_consumed':False,
      'arithmetic_boundary':'different integer cross-product code;Python arbitrary-precision integer runtime shared',
      'registration_sha256':sha(regpath),'code_sha256':sha(__file__),'source_manifest_sha256':sha(reg['input_manifest']['path'])}
    with out.open('x') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'pid':result['pid'],'checks':len(checks),'seconds':result['seconds']}))
if __name__=='__main__':main()
