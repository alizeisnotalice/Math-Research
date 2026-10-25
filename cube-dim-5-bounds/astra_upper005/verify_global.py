"""Consume an independently generated degree-26 space-moment artifact.

Input: first CLI argument, either the raw denominator/accumulator schema used by 004
or a rational space_moments.json with families d1,d2,b. No backend is run here.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import sys
import time

ROOT=Path(__file__).resolve().parent


def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    began=time.monotonic()
    if len(sys.argv)!=2:
        raise SystemExit('Provide exact backend moment JSON path; no implicit stale input is allowed.')
    backend=Path(sys.argv[1]).resolve()
    data=json.loads(backend.read_text())
    prior=load_module('upper004_global_locked',ROOT.parent/'astra_upper004'/'verify_global.py')
    helper=load_module('upper005_checkpoint_helpers',ROOT/'checkpoint.py')
    if 'families' in data:
        families={name:[{tuple(map(int,k.split(','))):F(v) for k,v in poly.items()} for poly in data['families'][name]] for name in ['d1','d2','b']}
    else:
        denominators=list(map(int,data['denominators']))
        assert len(denominators)>=27
        for n in range(1,27):assert denominators[n]==denominators[n-1]*(4*n-1)
        families={name:prior.read_family(data[name],denominators) for name in ['d1','d2','b']}
    assert len(families['d1'])>=15 and len(families['d2'])>=15 and len(families['b'])>=23
    oldspace_path=ROOT.parent/'astra_upper004'/'space_moments.json'
    oldspace=json.loads(oldspace_path.read_text())['families']
    for name in ['d1','d2','b']:
        for n,poly in enumerate(oldspace[name]):
            assert families[name][n]=={tuple(map(int,k.split(','))):F(v) for k,v in poly.items()}
    checkpoint_path=ROOT/'checkpoint.json'
    checkpoint=json.loads(checkpoint_path.read_text())
    keys={'d1':'E_d_W_power','d2':'E_d_squared_W_power','b':'E_b_corner_W_power'}
    for name,limit in [('d1',14),('d2',14),('b',22)]:
        for n,poly in enumerate(families[name][:limit+1]):
            assert all(sum(k)==n+(name=='b') for k in poly)
            assert prior.corner_value(poly)==F(checkpoint['moments_to_order_22'][keys[name]][n])
    portable=ROOT/'space_moments.json'
    if backend!=portable.resolve():
        portable_data={'families':{name:[{','.join(map(str,k)):str(v) for k,v in sorted(poly.items())} for poly in family] for name,family in families.items()},
                       'backend_source_sha256':data.get('backend_source_sha256'),
                       'raw_backend_sha256':hashlib.sha256(backend.read_bytes()).hexdigest(),
                       'canonical_moment_count':data.get('canonical_moment_count'),
                       'max_B_bits':data.get('max_B_bits')}
        portable.write_text(json.dumps(portable_data,indent=2)+'\n')
    A,B,c=(F(checkpoint['constants'][key]) for key in ['A_lower','B_upper','c_lower'])
    results=[]
    for candidate in checkpoint['quadratic_candidates']:
        if not candidate['passes_necessary_screen']:continue
        f=list(map(F,candidate['f_coefficients']))
        a=F(candidate['a_selected'])
        f7,f11=helper.ppow(f,7),helper.ppow(f,11)
        power={}
        for n,coef in enumerate(f7):
            factor=F(11,4)*a**7*coef
            for key,val in families['d1'][n].items():power[key]=power.get(key,F(0))+factor*A*val
            for key,val in families['d2'][n].items():power[key]=power.get(key,F(0))-factor*B*val
        for n,coef in enumerate(f11):
            factor=-F(7,4)*a**11*coef
            for key,val in families['b'][n].items():power[key]=power.get(key,F(0))+factor*val
        power={k:v for k,v in power.items() if v}
        assert prior.corner_value(power)==F(candidate['corner_value'])
        face,degrees,tensor,index=prior.bernstein_face(power)
        minimum=tensor[index]
        lower=c*minimum
        simple=F(lower.numerator*10**8//lower.denominator,10**8)
        results.append({'kappa':candidate['kappa'],'a':str(a),'f_coefficients':candidate['f_coefficients'],
                        'symmetric_power_coefficients':{','.join(map(str,k)):str(v) for k,v in sorted(power.items())},
                        'face_Bernstein_degrees':degrees,'coefficient_count':len(tensor),
                        'minimum_Bernstein_index':index,'minimum_Bernstein_value':str(minimum),
                        'positive_global_certificate':minimum>0,
                        'global_g_exact_rational':str(lower) if minimum>0 else None,
                        'global_g_simple_rational':str(simple) if minimum>0 else None,
                        'strictly_improves_7_over_156250':minimum>0 and simple>F(7,156250),
                        'float_display_only':float(lower),
                        'Bernstein_coefficients_sha256':hashlib.sha256(json.dumps({','.join(map(str,k)):str(v) for k,v in sorted(tensor.items())},sort_keys=True).encode()).hexdigest()})
    out={'status':'full_space_verification_of_three_registered_quadratic_candidates','new_C5_upper_bound':False,
         'backend_artifact':str(backend),'backend_sha256':hashlib.sha256(backend.read_bytes()).hexdigest(),
         'portable_space_moments_sha256':hashlib.sha256(portable.read_bytes()).hexdigest(),
         'independent_backend_source_sha256':data.get('backend_source_sha256'),
         'crosschecks':'all 004 space coefficients and all available independent 005 corner moments agree exactly',
         'candidates':results,
         'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'registration.md',ROOT/'checkpoint.py',checkpoint_path,Path(__file__),ROOT/'graph_moments.cpp',oldspace_path,ROOT.parent/'astra_upper004'/'verify_global.py']},
         'runtime':{'python':sys.executable,'seconds':time.monotonic()-began}}
    (ROOT/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'candidates':[{k:v for k,v in q.items() if k!='symmetric_power_coefficients'} for q in results],'runtime':out['runtime']},indent=2))


if __name__=='__main__':main()
