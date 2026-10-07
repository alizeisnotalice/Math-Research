"""Read saved outputs; re-evaluate only individual root physical scales."""
from pathlib import Path
import importlib.util,json,math
import numpy as np
D=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('model',D/'model.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
points={z['id']:z for z in json.loads((D/'frozen_points.json').read_text())['points']};records=[];cross=[]
for rd in [1,2,3]:
    z=json.loads((D/f'round{rd}.json').read_text())
    for row in z['rows']:
        point=points[row['point_id']];root=row['finite_visible_FIRST'];sig=root['sigma_raw'];L=root['Ls'];n=row['n'];N=row['spectral_N']
        h,b,e,c=m.coordinates(np.array(point['Z']),np.array(point['dx']),[L],N)
        data=np.load(D/f'round{rd}_{point["id"]}.npz');S=data['sigma'];G=data['response_mid'].max(axis=1);j=int(np.where(G>=1)[0][-1]);a=float(S[j]);d=float(S[j+1])
        endpoints=[];bounds=[]
        for bb in [np.maximum(0,b-e),b+e]:
            fa=float(m.response(h,bb,a)[0]);fd=float(m.response(h,bb,d)[0]);endpoints.append([fa,fd])
            if fa>1 and fd<1:
                lo=a;hi=d
                for _ in range(38):
                    mid=(lo+hi)/2
                    if m.response(h,bb,mid)[0]>=1:lo=mid
                    else:hi=mid
                bounds.append((lo+hi)/2)
            else:bounds.append(None)
        records.append(dict(round=rd,point_id=point['id'],physical_L=L,finite_visible_raw_sigma=sig,fixed_L_scalar_integral_uncertainty_root_lower=bounds[0],fixed_L_scalar_integral_uncertainty_root_upper=bounds[1],width=None if None in bounds else bounds[1]-bounds[0],endpoint_response_checks=endpoints,scope='Fixed finite-L candidate root only; analytic scalar quadrature interval, plus unverified rounding guard. Bisection digits do not control kernel quadrature or hidden continuous future crossings.'))
        if rd==3:
            alpha=math.ceil(math.sqrt(n));order=math.ceil((n+alpha)/2)
            nodes,weights=np.polynomial.legendre.leggauss(order);v=(nodes+1)/2;w=weights/2
            vals=np.array([m.response(h,b,sig+(1-sig)*vv)[0] for vv in v])
            gl=float(np.sum(w*alpha*v**(alpha-1)*vals));dp=m.weighted_future_average(h[:,0,:],b[:,0,:],sig)
            assert abs(dp-gl)<1e-11*max(1,abs(dp)),(point['id'],dp,gl)
            traces=data['hard_jump_traces'];ti=int(np.argmin(abs(traces[:,0]-L)))
            cross.append(dict(point_id=point['id'],n=n,alpha=alpha,Gauss_order=order,polynomial_degree=n+alpha-1,positive_Bernstein=dp,independent_polynomial_exact_Gauss=gl,abs_difference=abs(dp-gl),root_L_jump_trace=traces[ti].tolist(),root_L_jump_size=float(traces[ti,2]-traces[ti,1]),every_arrival_present_in_saved_grid=all(np.isin(point['L_breakpoints'],data['L']))))
(D/'root_error_brackets.json').write_text(json.dumps(dict(status='NUMERICAL_FIXED_L_QUADRATURE_UNCERTAINTY_BRACKETS',records=records),indent=2))
(D/'future_average_crosscheck.json').write_text(json.dumps(dict(status='PASS_POSITIVE_BERNSTEIN_VS_POLYNOMIAL_EXACT_GAUSS_FLOAT',records=cross,scope='Gauss exactness is a polynomial identity in exact arithmetic; numerical nodes/weights/function values remain floating, not outward certified.'),indent=2))
print('POSTCHECK PASS',max(z['abs_difference'] for z in cross))
