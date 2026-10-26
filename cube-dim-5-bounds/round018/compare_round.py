"""New exact round018-vs017 comparison; no roots or physical evaluations."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,os,sys,time
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    reg,out=map(Path,sys.argv[1:]);assert not out.exists();r=json.loads(reg.read_bytes());t=time.monotonic()
    for row in r['sources'].values():assert sha(row['path'])==row['sha256']
    old=json.loads(Path(r['sources']['round017']['path']).read_bytes())
    new=json.loads(Path(r['sources']['round018']['path']).read_bytes())
    oldI,newI=Q(old['complete_I_upper']),Q(new['complete_I_upper'])
    oldg,newg=Q(old['qualified_g']),Q(new['qualified_g'])
    old4,new4=Q(old['best_upper']['fourth_power']),Q(new['best_upper']['fourth_power'])
    ol,ou=map(Q,old['best_upper']['interval']);nl,nu=map(Q,new['best_upper']['interval'])
    assert 0<newI<oldI and newg>oldg>0 and 0<new4<old4<32**4
    assert nu<ol and new['lower_changed'] is False and new['dual_breakthrough_at'] is None
    result={'status':'PASS_complete_exact_round018_comparison','pid':os.getpid(),'seconds':time.monotonic()-t,
            'exact_I_decrease':str(oldI-newI),'exact_g_increase':str(newg-oldg),
            'exact_fourth_power_decrease':str(old4-new4),'upper_decrease_interval':[str(ol-nu),str(ou-nl)],
            'simple_upper_decrease':str(Q(old['simple_rational_upper'])-Q(new['simple_rational_upper'])),
            'comparisons_vs32':new['comparisons'],'L_ref_decrease':'0','lower_original_sha256':r['sources']['lower_original']['sha256'],
            'dual_breakthrough_at':None,'post_dual_rounds':0,'new_H_log_R_power_roots':0,'registration_sha256':sha(reg)}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'pid':result['pid'],'seconds':result['seconds']}))
if __name__=='__main__':main()
