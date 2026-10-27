"""One new integer-only audit of outward paper constants; no Fraction/root."""
from pathlib import Path
import hashlib,json,os,time
H=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_bytes())
def rational(s):
 p=s.split('/');assert len(p) in (1,2);n=int(p[0]);d=1 if len(p)==1 else int(p[1]);assert d>0;return n,d

def main():
 target=H/'result.json';assert not target.exists();reg=read(H/'registration.json');began=time.monotonic()
 for p,r in reg['sources'].items():assert Path(p).stat().st_size==r['bytes'] and sha(p)==r['sha256']
 a=read(reg['complete_certificate']);paper=read(reg['paper_constants']);N=int(paper['denominator']);G=int(paper['g_lower_numerator']);J=int(paper['I_upper_numerator']);K=28936079
 assert (G,J,N)==(11231606208775453,4474860503582229,10**20)
 gn,gd=rational(a['qualified_g']);inn,ind=rational(a['complete_I_upper']);cn,cd=rational(a['simple_rational_upper'])
 assert (cn,cd)==rational(paper['C5_upper'])==(K,10**6)
 assert rational(paper['g_lower'])==(G,N) and rational(paper['I_upper'])==(J,N)
 lower_residual=gn*N-G*gd;upper_residual=J*ind-inn*N
 assert G>0 and 0<=lower_residual<gd and 0<=upper_residual<ind
 assert G*gd<=gn*N and inn*N<=J*ind
 gap=K**4*G**5-3125*J**4*10**44;assert gap>0
 pn,pd=rational(paper['strict_polynomial_gap']);den=10**24*N**5
 assert pn>0 and pn*den==gap*pd
 assert paper['source']['path']==reg['complete_certificate'] and paper['source']['sha256']==reg['sources'][reg['complete_certificate']]['sha256']
 for p,r in reg['sources'].items():assert sha(p)==r['sha256']
 r={'status':'independent_integer_outward_paper_scalar_audit_pass','pid':os.getpid(),'seconds':time.monotonic()-began,'registration_sha256':sha(H/'registration.json'),'G':G,'H':J,'N':N,'K':K,'g_floor_crossmultiply_residual':str(lower_residual),'I_ceil_crossmultiply_residual':str(upper_residual),'strict_integer_gap':str(gap),'gap_identity_denominator':str(den),'from_original_exact_certificate':True,'same_existing_portable_C5_bound':True,'new_root_H_R_mean_log_calls':0,'old_successful_script_replays':0,'frozen178_report_modified':False,'global_G180_status_unchanged':True}
 with target.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({'status':r['status'],'pid':r['pid'],'seconds':r['seconds'],'result_sha256':sha(target)}))
if __name__=='__main__':main()
