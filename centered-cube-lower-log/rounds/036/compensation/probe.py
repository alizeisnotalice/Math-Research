"""Frozen random-source Fourier projection and legal-density budgets.

Finite pair formulas include exact box-smoothing algebra; all evaluations
are ordinary NumPy floating point, not interval certificates.
"""
from pathlib import Path
import hashlib,json,math,sys
import numpy as np

HERE=Path(__file__).resolve().parent
PLAN={'seed':202610290361,'K':[63,255,1023],'v':[8,32,128],
      'p':.5,'w':4,'L_grid':np.linspace(10/9,11/9,129).tolist(),
      'comparison':'original, exact-low-projection, H-rescaled, B-calibrated',
      'centering':'Remove primitive spatial mean before measured energy normalization.',
      'scope':'Source diagonal/positivity diagnostics only; no MC tails or C_n claim.'}
def save(name,value):
    (HERE/name).write_text(json.dumps(value,indent=2)+'\n')
def pair_data(a):
    d=(a[:,None]-a[None,:])%1
    d=np.sort(d[d>0])
    return d,np.r_[0.,np.cumsum(d)]
def smooth_diagonal(k,u,t,delta,d,cs):
    ix=np.searchsorted(d,t,side='right')
    un=k*t/2-k*k*t*t/2+t*ix-cs[ix]
    correction=[]
    for tt in t:
        near=d[(d>tt-delta)&(d<tt+delta)]
        x=tt-near
        hinge=((x+delta)**3-2*np.maximum(x,0)**3)/(6*delta*delta)
        correction.append(float(np.sum(hinge-np.maximum(x,0))))
    short=d[d<delta]
    short_extra=float(np.sum((delta-short)**3)/(3*delta*delta))
    expected_DV=k*delta/6-k*k*delta*delta/12+short_extra
    r0=float(np.mean(u*u-u+1/3))
    dsm=un-k*delta/6+np.array(correction)-short_extra
    return 2*dsm,r0-expected_DV,un*2

def main():
    if (HERE/'results.json').exists()or(HERE/'registration.json').exists():
        raise RuntimeError('Refusing repeated experiment.')
    plan=dict(PLAN,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    save('registration.json',plan)
    rng=np.random.default_rng(PLAN['seed']);L=np.array(PLAN['L_grid']);t=L%1
    rows=[];sources=[]
    for k in PLAN['K']:
        u=rng.random(k);a=(np.arange(k)+u)/k
        sources.append({'K':k,'U':u.tolist()})
        d,cs=pair_data(a);meanp=float(.5-u.mean())
        for v in PLAN['v']:
            M=math.ceil(math.sqrt(v));m=np.arange(1,M+1);delta=1/(4*k*v)
            b=4*(729/29)*((np.sin(np.pi*m*10/27)/(np.pi*m))**4
                         +3*(np.sin(np.pi*m/9)/(np.pi*m))**4)
            c=np.sqrt(2*b/3);gamma=np.pi*m*c
            z=np.exp(-2j*np.pi*a[:,None]*m).sum(axis=0)
            zd=z*np.sinc(m*delta);ph=zd/(2j*np.pi*m)
            hh,r0sm,hh_un=smooth_diagonal(k,u,t,delta,d,cs)
            s2=np.sin(np.pi*L[:,None]*m)**2
            bb=(s2*(c*c/2)).sum(axis=1)
            cross=(s2*(-4*c*ph.imag)).sum(axis=1)
            removed=(s2*(8*abs(ph)**2)).sum(axis=1)
            hh_high=hh-removed
            C0=r0sm-meanp**2-2*float(np.sum(abs(ph)**2))
            assert C0>0
            hscale=1/math.sqrt(6*C0);bscale=math.sqrt(6*C0)
            target=9*L*L/29
            gram={'original':hh+bb+cross,
                  'projected':hh_high+bb,
                  'H_rescaled':hscale*hscale*hh_high+bb,
                  'B_calibrated':hh_high+(6*C0)*bb}
            targets={name:(target*(6*C0)if name=='B_calibrated'else target)for name in gram}
            p=.5
            lower={
                'original':1-p-p/k*float(np.sum(abs(gamma))),
                'projected':1-p-p/k*float(np.sum(abs(gamma-2*zd))),
                'H_rescaled':1-p*hscale-p/k*float(np.sum(abs(gamma-2*hscale*zd))),
                'B_calibrated':1-p-p/k*float(np.sum(abs(bscale*gamma-2*zd)))}
            row={'K':k,'v_reference':v,'M':M,'delta':delta,
                 'primitive_mean':meanp,'R0_smoothed':r0sm,'C0_centered_high':C0,
                 'H_scale':hscale,'B_scale':bscale,'target_scale_B_calibrated':6*C0,
                 'effective_v_B_calibrated':6*C0*v,
                 'peak_mass_H_rescaled':p*hscale,
                 'density_triangle_lower':lower,'all_lower_bounds_positive':min(lower.values())>0,
                 'density_mass_exact_by_formula':1,
                 'projection_derivative_envelope':float(2*np.sum(abs(zd))),
                 'P_low_sup_envelope':float(2*np.sum(abs(ph))),
                 'H_high_sup_envelope':2+4*float(np.sum(abs(ph))),
                 'HB_original_max_abs':float(abs(cross).max()),
                 'HB_projected_exact':0,'removed_HH_max':float(removed.max()),
                 'smoothing_HH_max_difference':float(abs(hh-hh_un).max()),
                 'metrics':{name:{'max_abs_Gram_error':float(abs(g-targets[name]).max()),
                    'v_times_max_error':float(v*abs(g-targets[name]).max()),
                    'max_relative_Gram_error':float((abs(g-targets[name])/targets[name]).max()),
                    'min_Gram':float(g.min()),
                    'diagonal_error_spread':float(np.ptp(g-targets[name]))}
                    for name,g in gram.items()},
                 'curves':{'L':L.tolist(),'HH_original':hh.tolist(),'HH_high':hh_high.tolist(),
                           'BB':bb.tolist(),'cross_original':cross.tolist(),
                           'Gram':{name:g.tolist()for name,g in gram.items()},
                           'target':{name:g.tolist()for name,g in targets.items()}}}
            rows.append(row)
    save('source_jitters.json',sources)
    save('results.json',{'registration_SHA256':hashlib.sha256((HERE/'registration.json').read_bytes()).hexdigest(),
        'python':sys.executable,'numpy':np.__version__,'rows':rows,
        'caution':'Ordinary float finite-grid source diagnostics; positivity triangle margins not interval certified; exact orthogonality and mass refer to the analytic definitions; no tail or C_n conclusion.'})
    for row in rows:
        print(json.dumps({'K':row['K'],'v':row['v_reference'],
            'density_lb':min(row['density_triangle_lower'].values()),
            'v_error':{n:round(x['v_times_max_error'],6)for n,x in row['metrics'].items()}},sort_keys=True))

if __name__=='__main__':main()
