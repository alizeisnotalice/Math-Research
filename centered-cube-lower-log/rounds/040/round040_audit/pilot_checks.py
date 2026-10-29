"""Dynamic-R interval checks and independent exact tiny count-proposal/RN checks."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,sys,json,math,hashlib,itertools
import numpy as np
HERE=Path(__file__).resolve().parent;OTHER=HERE.parent/'round040_stream'
sys.dont_write_bytecode=True;sys.path.insert(0,str(OTHER))
spec=importlib.util.spec_from_file_location('stream40',OTHER/'probe.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
dynamic=[]
for R in [4,8,16,32]:
    src=mod.Source(R);cuts,pl,pu=mod.make_cuts(src);J=3;N=2
    mu=np.full((N,R,J),src.n*src.pchi)
    var=np.full_like(mu,src.n*src.pchi*(1-src.pchi))
    z=np.zeros_like(mu);lr=np.zeros_like(mu)
    err=np.full((R,J),src.n*1e-13)
    bounds,diag=mod.bounds(np.array([z,mu,var,lr]),src,cuts,pl,pu,np.zeros((R,J)),err)
    assert all(np.array_equal(a,b)for a,b in zip(bounds['joint'],bounds['prior']))
    assert len(cuts)==R-1 and len(set(cuts))==R-1
    assert pl.sum()<=1<=pu.sum() and pl.min()>0
    dynamic.append({'R':R,'cuts_count':len(cuts),'prior_lower_sum':float(pl.sum()),
                    'prior_upper_sum':float(pu.sum()),'minimum_prior_lower':float(pl.min()),
                    'identical_PB_and_prior_mock_intervals_match':True})
# Three-state original query count law and four normalized product tilts.
n=3;Q=3;R=2;J=2;components=[list(map(F,x))for x in [
    ['3/5','1','7/5'],['4/5','7/5','8/5'],['1/2','4/5','7/5'],['6/5','9/10','6/5']]]
norm=[sum(a)/Q for a in components]
CE=F(101,100);CZ=max([F(1)]+[a**n for a in norm]);cap=2*R*J*CZ*CE
totalQ=F(0);weighted=F(0);checks=0
for c0 in range(n+1):
 for c1 in range(n-c0+1):
    counts=[c0,c1,n-c0-c1];multi=math.factorial(n)//math.prod(math.factorial(c)for c in counts)
    P0=F(multi,Q**n)
    M=[math.prod(a**c for a,c in zip(comp,counts))for comp in components]
    D=F(1,2)+sum((m/z**n for m,z in zip(M,norm)),F(0))/(2*R*J)
    qprob=P0/2+sum((multi*math.prod((a/(Q*z))**c for a,c in zip(comp,counts))
                       for comp,z in zip(components,norm)),F(0))/(2*R*J)
    assert qprob/P0==D and 1/D<=2
    totalQ+=qprob;weighted+=qprob/D
    for lam in [F(1,2),F(2),F(5)]:
        # Markov/union envelopes for both joint and prior mean upper arrays.
        U=min(F(1),CE*sum(M)/lam)
        assert lam*U/D<=min(2*lam,cap)
        B=U/2;assert lam*B/D<=min(2*lam,cap);checks+=2
assert totalQ==weighted==1
out={'dynamic_R_checks':dynamic,'exact_count_RN':{'n':n,'Q':Q,'R':R,'J':J,
        'proposal_mass':str(totalQ),'E_proposal_inverse_RN':str(weighted),
        'two_direction_contribution_cap_checks':checks,'normalizers':list(map(str,norm)),
        'C_Z':str(CZ),'C_E':str(CE)},
     'CI_family':2*len(mod.RVALUES)*len(mod.MODELS)*len(mod.H),
     'alpha_each':mod.ALPHA,'no_MC_or_large_cache':True,
     'reviewed_SHA256':{name:hashlib.sha256((OTHER/name).read_bytes()).hexdigest()for name in ['probe.py','helpers.py']}}
assert out['CI_family']==504
target=HERE/'pilot_checks.json'
if target.exists():raise RuntimeError('Refusing replacement')
target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
