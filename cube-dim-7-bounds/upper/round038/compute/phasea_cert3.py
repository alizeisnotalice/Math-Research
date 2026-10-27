#!/usr/bin/env python3
"""At most three selected Phase A patterns: exact primal/dual feasibility certificates."""
import os,sys,json,time,hashlib,resource
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
from joint_core import rational_bounds,maskids,load
P=Path(__file__).resolve();D=P.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def select():
 reg=json.loads((D/'phasea_registration.json').read_text());summary=json.loads((D/'phasea_summary.json').read_text());assert len(summary['results'])==79
 order={(r['template'],r['stratum'],r['caseindex'],r['pattern']):i for i,r in enumerate(reg['identities'])};rows=[]
 for rec in summary['results']:
  path=D/rec['output'];assert sha(path)==rec['sha256'];r=json.loads(path.read_text());id=r['historical_identity'];idx=order[(id['template'],id['stratum'],id['caseindex'],id['pattern'])];L=Q(r['L_after']);L0=Q(r['L_before']);Khi=Q(r['U'])/Q(r['history_M_lower_unverified']);rows.append({'path':path,'index':idx,'stop':r['stop'],'residual':(L-Khi)/L,'reduction':(L0-L)/L0})
 chosen=[]
 def pick(candidates,key):
  if candidates:
   v=max(candidates,key=lambda r:(r[key],-r['index']));chosen.append(v)
 pick([r for r in rows if r['stop']=='stop_pair_coordinate' and r['residual']>0],'residual')
 pick([r for r in rows if r not in chosen and r['reduction']>0],'reduction')
 pick([r for r in rows if r not in chosen and r['residual']>0],'residual')
 return chosen

def register():
 chosen=select();out={'batch':'R038 Phase A independent LP certification','rules':'SPEC §7: max stationary historical residual, max relative reduction, max remaining historical residual; manifest-order ties','limit':3,'source_sha256':sha(P),'core_sha256':sha(D/'joint_core.py'),'phasea_summary_sha256':sha(D/'phasea_summary.json'),'selected':[{'file':r['path'].name,'sha256':sha(r['path']),'manifest_index':r['index'],'historical_residual':str(r['residual']),'relative_reduction':str(r['reduction'])} for r in chosen]};(D/'phasea_cert3_registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

def run():
 reg=json.loads((D/'phasea_cert3_registration.json').read_text());assert reg['source_sha256']==sha(P) and reg['core_sha256']==sha(D/'joint_core.py') and reg['phasea_summary_sha256']==sha(D/'phasea_summary.json');outputs=[]
 for sel in reg['selected']:
  t=time.monotonic();cpu=time.process_time();path=D/sel['file'];assert sha(path)==sel['sha256'];r=json.loads(path.read_text());h=list(map(Q,r['witness_budgets']));C=[sum((1<<i for i in c),0) for c in r['maximal_cliques']];U=Q(r['U']);q=list(map(Q,r['q_after']));L0=Q(r['L_before']);Lf=load(q,C);assert Lf==Q(r['L_after']) and sum((h[i]*q[i] for i in range(len(h))),Q())==U
  Mlo,Mhi,z,y=rational_bounds(h,C);assert sum(z,Q())==Mhi and sum((h[i]*y[i] for i in range(len(h))),Q())==Mlo
  assert all(sum((z[j] for j,c in enumerate(C) if c>>i&1),Q())>=h[i] for i in range(len(h))) and all(sum((y[i] for i in maskids(c)),Q())<=1 for c in C)
  klo=U/Mhi;khi=U/Mlo;eligible=L0-khi>=L0/10000
  verdict='not_certified_eligible' if not eligible else ('H_pass' if 2*Lf<=L0+klo else 'H_fail' if 2*Lf>L0+khi else 'H_unresolved')
  out={'source':sel,'M_lower':str(Mlo),'M_upper':str(Mhi),'primal_m':[str(v) for v in z],'dual_y':[str(v) for v in y],'primal_feasible':True,'dual_feasible':True,'initial_gap_certified_eligible':eligible,'Lopt_lower':str(klo),'Lopt_upper':str(khi),'L0':str(L0),'Lfinal':str(Lf),'verdict':verdict,'pid':os.getpid(),'wall_seconds':time.monotonic()-t,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
  dest=D/('lp_cert_'+sel['file']);dest.write_text(json.dumps(out,indent=2)+'\n');outputs.append({'file':dest.name,'sha256':sha(dest),'verdict':verdict})
 (D/'phasea_cert3_results.json').write_text(json.dumps({'registration_sha256':sha(D/'phasea_cert3_registration.json'),'outputs':outputs},indent=2)+'\n');print(json.dumps(outputs,indent=2))
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register()
 elif sys.argv[1:]==['--run']:run()
 else:raise ValueError('use --register or --run')
