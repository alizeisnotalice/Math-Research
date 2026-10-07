#!/usr/bin/env python3
"""Original c=1 diagonal Levy quadrature diagnostic; no old imports/reruns."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from numpy.polynomial.legendre import leggauss

BASE = Path(__file__).resolve().parent
PREFIX = "diagonal_levy_spatial_probe_20261007"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gauss_unit(nodes):
    t,w = leggauss(nodes)
    return (t+1)/2,w/2


def spectral_rho(a,s):
    ell = np.log(s*s-1)
    denominator = (a*s*s-(1-a)*ell)**2 + np.pi**2*(1-a)**2
    return a*s*s/denominator


def phase(r,s):
    u = (1-r)/r*s*s-np.log(s*s-1)
    return np.arctan2(np.pi,u)/(2*np.pi)


def quadrature(s_nodes):
    t,w=gauss_unit(s_nodes)
    s=1+t/(1-t)
    return s,w/(1-t)**2


def levy_profile(r,L,x,s,w):
    return np.exp(-np.outer(x/L,s)) @ (w*phase(r,s))/L


def generator_profile(r,q,x,s,w):
    exponent=np.exp(-np.outer(x,s))
    f=phase(r,s)
    # q dr Lambda - Lambda - x Lambda': evaluate nonsingular original formula.
    h=q*spectral_rho(1-r,s)/(2*(1-r))
    return exponent @ (w*(h-f)) + x*(exponent @ (w*s*f))


def advance(r,kind,step):
    if kind == 'root':
        return 1-(1-r)*np.exp(-2*step)
    factor=float(kind)
    odds=r/(1-r)*np.exp(2*factor*step)
    return odds/(1+odds)


def qvalue(r,kind):
    return 2*(1-r) if kind=='root' else float(kind)*2*r*(1-r)


def main():
    regpath=BASE/(PREFIX+'_registration.json')
    reg=json.loads(regpath.read_text())
    resultpath=BASE/(PREFIX+'_results.json')
    assert not resultpath.exists(), 'Do not overwrite old results.'
    start=time.perf_counter()
    rounds=[]
    rows=[]
    validations=[]
    profile_files=[]
    for plan in reg['rounds']:
        s,w=quadrature(plan['s_gauss_nodes'])
        x=np.r_[0,np.geomspace(1e-6,1,plan['x_nodes']-1)]
        aq,aw=gauss_unit(plan['a_gauss_nodes'])
        profiles={'x':x,'s':s,'s_weights':w}
        indices=[]
        for r in reg['r_values']:
            a=1-r
            rho=spectral_rho(a,s)
            for xi in [0,.5,2,8]:
                approx=np.sum(w*rho*2*s/(s*s+xi*xi))
                if xi==0:
                    exact=1.
                else:
                    v=xi*xi
                    exact=np.log1p(v)/(a*v+(1-a)*np.log1p(v))
                validations.append(dict(round=plan['name'],contract='original_Ga_symbol',r=r,a=a,xi=xi,
                    quadrature_value=float(approx),original_symbol=float(exact),residual=float(approx-exact)))
            aa=a+(1-a)*aq
            for sv in [1.0001,1.1,2.,8.]:
                approx=.5*(1-a)*np.sum(aw*spectral_rho(aa,sv)/aa)
                closed=float(phase(r,sv))
                validations.append(dict(round=plan['name'],contract='a_integral_closed_phase',r=r,s=sv,
                    quadrature_value=float(approx),closed_value=closed,residual=float(approx-closed)))
            old=levy_profile(r,1,x,s,w)
            for kind in [.25,.5,1.,2.,'root']:
                q=qvalue(r,kind)
                g=generator_profile(r,q,x,s,w)
                key=f'r{r}_slope{kind}'
                profiles[key+'_generator']=g
                record=dict(round=plan['name'],r=r,slope_kind=str(kind),q=float(q),
                    minimum_generator=float(g.min()),minimum_x=float(x[np.argmin(g)]),
                    negative_part_grid_integral=float(2*np.trapezoid(np.maximum(-g,0),x)),
                    analytic_positive_sufficient=kind=='root' or float(kind)>=1,
                    finite_steps=[])
                for step in reg['logL_steps']:
                    newr=advance(r,kind,step)
                    delta=levy_profile(newr,float(np.exp(step)),x,s,w)-old
                    profiles[key+f'_step{step}']=delta
                    record['finite_steps'].append(dict(logL_step=step,new_r=float(newr),
                        minimum_difference=float(delta.min()),minimum_x=float(x[np.argmin(delta)]),
                        negative_part_grid_integral=float(2*np.trapezoid(np.maximum(-delta,0),x)),
                        total_mass_difference_exact=.5*float(np.log((1-r)/(1-newr)))))
                indices.append(len(rows)); rows.append(record)
        profpath=BASE/(PREFIX+'_'+plan['name']+'_profiles.npz')
        assert not profpath.exists()
        np.savez_compressed(profpath,**profiles)
        profile_files.append(dict(path=str(profpath),sha256=sha(profpath)))
        rounds.append(dict(**plan,row_indices=indices))
    comparisons=[]
    for i in range(2):
        left,right=rounds[i],rounds[i+1]
        for li,ri in zip(left['row_indices'],right['row_indices']):
            l,r=rows[li],rows[ri]
            comparisons.append(dict(coarse=left['name'],fine=right['name'],r=r['r'],slope_kind=r['slope_kind'],
                minimum_generator_change=r['minimum_generator']-l['minimum_generator'],
                negative_integral_change=r['negative_part_grid_integral']-l['negative_part_grid_integral']))
    payload=dict(status='complete',registration_sha256=sha(regpath),script_sha256=sha(Path(__file__)),
        rounds=rounds,records=rows,validations=validations,refinement_diagnostics=comparisons,
        profile_files=profile_files,elapsed_seconds=time.perf_counter()-start,
        no_clipping=True,no_confidence_interval=True,no_continuous_numerical_certificate=True,
        analytic_positive_contract='q>=2r(1-r); proved in companion note, not inferred from samples',
        scope='original c=1 one-coordinate kernel; tensor applicability analytic; no actual FIRST/geom pressure')
    resultpath.write_text(json.dumps(payload,indent=2)+'\n')
    summary={p['name']:dict(records=sum(row['round']==p['name'] for row in rows),
        positive_candidate_min=min(row['minimum_generator'] for row in rows if row['round']==p['name'] and row['analytic_positive_sufficient']),
        max_symbol_residual=max(abs(v['residual']) for v in validations if v['round']==p['name'] and v['contract']=='original_Ga_symbol')) for p in reg['rounds']}
    print(json.dumps(dict(elapsed_seconds=payload['elapsed_seconds'],summary=summary),indent=2))
    print(str(resultpath))


if __name__=='__main__':
    main()
