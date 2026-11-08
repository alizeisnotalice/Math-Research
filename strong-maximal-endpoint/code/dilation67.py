"""Even-box angular means: candidate extrema plus an analytical uniform bound."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import math,json,time,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import irwinhall
from scipy.optimize import minimize_scalar
import scipy

def one(sp):
    b,D=sp
    if D==0:return {'b':b,'D':D,'sum_candidate':.5,'log_angular_fraction_candidate':0.,'solver_status':'exact lower-slice plateau','solver_success':True}
    def neg(s):
        p=irwinhall.logpdf(s,b)
        return -(math.lgamma(b)+p-(b-1)*math.log(b*D+s))
    z=minimize_scalar(neg,bounds=(1e-10,b-1e-10),method='bounded',options={'xatol':1e-9})
    return {'b':b,'D':D,'sum_candidate':float(z.x),'log_angular_fraction_candidate':float(-z.fun),'solver_status':str(z.message),'solver_success':bool(z.success)}
if __name__=='__main__':
    cfg=json.loads(Path(sys.argv[1]).read_text())['dilation'];outdir=Path(sys.argv[2]);outdir.mkdir(parents=True,exist_ok=True);workers=int(sys.argv[3]) if len(sys.argv)>3 else 12
    st=time.time();specs=[(b,D) for b in cfg['block_sizes'] for D in cfg['D']]
    with ProcessPoolExecutor(max_workers=workers) as ex:base=list(ex.map(one,specs))
    idx={(r['b'],r['D']):r for r in base};rows=[]
    for n in cfg['N']:
        a=1+1/n**2
        for r in cfg['block_counts']:
            if n//r<2 or n%r:continue
            b=n//r
            for D in cfg['D']:
                z=idx[b,D];logB=math.log(a)+r*z['log_angular_fraction_candidate'];upper=math.log(a)-(n-r)*math.log(D) if D else None
                if upper is not None:assert logB<=upper+1e-8
                rows.append({'N':n,'blocks':r,'block_size':b,'D':D,'amplitude':a,'log_mean_candidate':logB,'analytical_log_upper':upper,'log_Jensen_modular_ratio':(n-1)*math.log1p(math.log1p(1/n**2)),'strict_input_margin':1/n**2,'radialized_halo_empty_certified_by_bound':D>=2,'log_witness_volume':n*(math.log(2*(1-math.exp(-1)))-D)})
    out={'status':'complete','scipy_version':scipy.__version__,'workers':workers,'python_version':sys.version,'seconds':time.time()-st,'base_cases':base,'rows':rows,'scope':'Extrema are floating candidates, not global numerical certificates; B_D<=a D^{-(N-r)} is a separate analytical bound. Even-box source; no new C_N bound.'}
    (outdir/'dilation67.json').write_text(json.dumps(out,indent=2));print('Completed dilation checks',len(rows),flush=True)
