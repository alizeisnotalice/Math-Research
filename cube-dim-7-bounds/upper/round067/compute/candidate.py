from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
p=Path(__file__).resolve().parent
def logbound(x):
 y=(x-1)/(x+1);assert 0<=y<1;s=2*sum((y**(2*k+1)/F(2*k+1) for k in range(48)),F(0));return s,s+2*y**97/(97*(1-y*y))
e=F(1,2**20);ea=F(1,65536);er=F(1,4096);rlo=1+F(1,2**30);rhi=F(2751,2000)-F(1,2**30);P0=1-e-F(1,2**60);ul=P0*(1-e)/(1+e);uh=(1+e)/(P0*(1-e));M=2*(1+e)*uh;b=2*(1-er)*(1-e)*ul-2*e;m=F(637,2000);gamma=(1-e)*m*b/M;Llo,Lhi=logbound(F(6,5));rL,_=logbound(rhi/rlo);mu0up=F(1,2)-F(5,2)*Llo;bridge=F(2211,50000)/(1-e);chi=F(9,200);sigma=F(83,1000);eta=F(83,1834);raw=F(317,1000)-F(83,917)
vals={'logmass_lower':(1-ea)*rL,'m':m,'gamma':gamma,'gamma_margin':gamma-F(317,1000),'mu0_upper':mu0up,'bridge_upper':bridge,'pin_gap':chi-F(223,5000)-F(1,5000),'eta_minus_chi':eta-chi,'tau':sigma/2-chi*(1-sigma),'raw':raw,'raw_minus_9_40':raw-F(9,40),'variance_upper':F(18233,100000)**3/(1-e)}
checks={'mass_gt_m':vals['logmass_lower']>m,'gamma_gt_clean':gamma>F(317,1000),'mu0_lt_221_5000':mu0up<F(221,5000),'bridge_lt_223_5000':bridge<F(223,5000),'pin_gap_identity':vals['pin_gap']==F(1,5000),'eta_identity':eta-chi==F(47,183400),'tau_identity':vals['tau']==F(47,200000),'raw_identity':raw==F(207689,917000) and raw-F(9,40)==F(341,229250),'variance_lt_1_150':vals['variance_upper']<F(1,150)}
out={'scope':'conditional scalar parameters; no target profile or operator theorem','method':'48-term atanh log series with upper remainder 2*y^97/(97*(1-y²))','code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'verdict':checks,'exact':{k:str(v) for k,v in vals.items()},'decimal':{k:float(v) for k,v in vals.items()},'executed_unix':time.time()};f=p/'candidate_certificate.json';f.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'verdict':checks,'decimal':out['decimal'],'source_sha256':out['code_sha256'],'certificate_sha256':hashlib.sha256(f.read_bytes()).hexdigest()},indent=2))
