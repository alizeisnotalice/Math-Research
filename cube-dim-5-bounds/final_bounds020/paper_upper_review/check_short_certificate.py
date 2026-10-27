"""One newly supplied paper display certificate, exact integer cross-products."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,sys,time
def read(p):return json.loads(Path(p).read_bytes())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ratio(s):
    z=s.split('/');assert len(z) in (1,2)
    n,d=int(z[0]),int(z[1]) if len(z)==2 else 1
    assert n>0 and d>0;return n,d
def main():
    rp,out=map(Path,sys.argv[1:]);assert not out.exists();reg=read(rp)
    assert reg['single_new_display_case']==1 and reg['code_sha256']==sha(__file__)
    before={k:sha(v['path']) for k,v in reg['sources'].items()}
    assert all(before[k]==v['sha256'] for k,v in reg['sources'].items())
    begin=time.monotonic();paper=read(reg['sources']['paper_constants']['path']);old=read(reg['sources']['accepted175']['path'])
    G,H,N,K=(paper[k] for k in ['G','H','N','K'])
    assert all(type(v) is int for v in [G,H,N,K])
    assert (G,H,N,K)==(11231606208775453,4472478819809024,10**20,28920678)
    assert paper['source_sha256']==before['accepted175']
    gn,gd=ratio(old['qualified_g']);In,Id=ratio(old['complete_I_upper'])
    rg=gn*N-G*gd
    ri=H*Id-In*N
    gap=K**4*G**5-3125*H**4*10**44
    assert rg>=0 and ri>=0 and gap>0
    assert str(gap)==paper['strict_integer_gap']
    after={k:sha(v['path']) for k,v in reg['sources'].items()};assert before==after
    result={'status':'PASS_one_new175_display_certificate','pid':os.getpid(),
      'utc':datetime.now(timezone.utc).isoformat(),'seconds':time.monotonic()-begin,
      'G':G,'H':H,'N':N,'K':K,'checks':{'G_over_N_le_original_g':True,
      'H_over_N_ge_original_complete_I175':True,'strict_integer_gate_positive':True},
      'g_crossmultiply_residual':str(rg),'I_crossmultiply_residual':str(ri),
      'strict_integer_gap':str(gap),'saved_primary_gap_matches':True,'source_hashes_before_after_equal':True,
      'new_display_cases':1,'new_bounds_searches':0,'root_H_R_mean_slope_gain_DFS_or_old39check_replays':0,
      'research_still_paused':True,'registration_sha256':sha(rp),'code_sha256':sha(__file__)}
    with out.open('x') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'pid':result['pid'],'seconds':result['seconds'],'new_display_cases':1}))
if __name__=='__main__':main()
