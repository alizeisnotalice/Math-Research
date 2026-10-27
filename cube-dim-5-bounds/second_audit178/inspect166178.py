"""New second-audit exact field checks; no root isolation or old script import."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,os,time
H=Path(__file__).resolve().parent
D=Path('/Users/zhengzhihao/Desktop/T')
MAIN=D/'Math-Research-cube-dim-5/cube-dim-5-bounds'
AUDIT=D/'Math-Research-cube5-audit136-sol/cube-dim-5-bounds/audit_combine169'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def get(p):return json.loads(Path(p).read_bytes())
def main():
 out=H/'inspection166_result.json';assert not out.exists();began=time.monotonic();reg=get(H/'inspection166_registration.json')
 for p,r in reg['sources'].items():assert sha(p)==r['sha256'] and Path(p).stat().st_size==r['bytes']
 a=get(MAIN/'round019/primary_upper166_acceptance.json');old=get(MAIN/'round018/primary_integral156_acceptance.json');prod=get(D/'cube-finite-dim-scratch/upper_combine166/run001/completion.json');audit=get(AUDIT/'run001/results.json');gcert=get(MAIN/'round018/primary_curvature155_acceptance.json')
 assert audit['worker_exitcodes']==[0]*15 and audit['source_hashes_before_after_equal'] and audit['new_gain_cases']==97971 and audit['nodes']==97972 and audit['roots']==52360
 uh=Q(a['new_UH']);assert uh==Q(prod['new_UH'])==Q(audit['new_UH'])==Q(3994994121862399564161069833676547838746859,91343852333181432387730302044767688728495783936)
 g=Q(a['qualified_g']);J=Q(a['same_complete_J139_upper']);I=uh+J
 assert g==Q(old['qualified_g'])==Q(gcert['qualified_g']) and J==Q(old['same_complete_J139_upper']) and I==Q(a['complete_I_upper'])<Q(old['complete_I_upper'])
 fourth=3125*I**4/g**5;assert fourth==Q(a['best_upper']['fourth_power'])<Q(old['best_upper']['fourth_power'])<32**4
 lo,hi=map(Q,a['best_upper']['interval']);ol,ou=map(Q,old['best_upper']['interval']);assert 0<lo<hi<ol<ou and lo**4<=fourth<hi**4 and hi<=Q(a['simple_rational_upper'])==Q(28936079,1000000)<ol
 assert Q(a['exact_I_gain_vs156'])==Q(old['new_UH'])-uh and Q(a['exact_fourth_power_gain_vs156'])==Q(old['best_upper']['fourth_power'])-fourth
 m=Q(gcert['minimum_margin110']);terms=list(map(Q,gcert['Taylor_terms']));assert m>0 and len(terms)==9 and terms[0]==1
 for k in range(1,9):assert k*terms[k]==m*terms[k-1]
 p8=sum(terms,Q(0));assert g<=Q(gcert['g0'])*p8<g+Q(1,1<<192)
 assert gcert['H128_floor_domain']==a['H128_floor_domain']=='complete_F_only' and a['G_floor_domain']=='complete_Q_via_G_homogeneity'
 for p,r in reg['sources'].items():assert sha(p)==r['sha256']
 r={'status':'new_second_audit_exact_field_comparisons_pass','pid':os.getpid(),'seconds':time.monotonic()-began,'registration_sha256':sha(H/'inspection166_registration.json'),'UH166':str(uh),'I166':str(I),'g192':str(g),'J139':str(J),'upper_root_interval_checked_by_powers':[str(lo),str(hi)],'simple_upper':'28936079/1000000','strict_U166_below_U156':True,'saved_Taylor8_recurrences_and_down192_checked':True,'new_root_isolations':0,'new_H_DD_log_gradient_or_full_tree_calls':0,'old_successful_script_replays':0,'lower_changed':False}
 out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'status':r['status'],'pid':r['pid'],'seconds':r['seconds'],'result_sha256':sha(out)}))
if __name__=='__main__':main()
