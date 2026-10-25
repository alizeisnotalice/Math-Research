"""One fixed K3 candidate: exact full-face polynomial and Bernstein certificate."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib
import importlib.util
import json
import sys
import time

ROOT=Path(__file__).resolve().parent


def load(path):
    spec=importlib.util.spec_from_file_location('locked_004_bernstein',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def decode(families):
    return {name:[{tuple(map(int,k.split(','))):F(v) for k,v in poly.items()} for poly in family] for name,family in families.items()}


def main():
    began=time.monotonic()
    if len(sys.argv)!=2:raise SystemExit('Provide raw backend JSON or portable space_moments.json.')
    input_path=Path(sys.argv[1]).resolve();data=json.loads(input_path.read_text())
    helper_path=ROOT.parent/'astra_upper004'/'verify_global.py';helper=load(helper_path)
    names=['d1','d2','d3','d4','b']
    if 'families' in data:
        families=decode(data['families'])
    else:
        assert data['completed_W_power']==11
        denominators=list(map(int,data['denominators']))
        assert len(denominators)==28 and denominators[0]==1
        for n in range(1,28):assert denominators[n]==denominators[n-1]*(4*n-1)
        families={name:helper.read_family(data[name],denominators) for name in names}
    assert set(families)==set(names)
    for name in names:assert len(families[name])==(12 if name=='b' else 8)
    old_path=ROOT.parent/'astra_upper004'/'space_moments.json'
    old=decode(json.loads(old_path.read_text())['families'])
    for name in ['d1','d2','b']:assert families[name]==old[name]
    corner_path=ROOT.parent/'astra_upper006'/'verification.json'
    corner=json.loads(corner_path.read_text())
    for k in range(1,5):
        for n,poly in enumerate(families['d'+str(k)]):
            assert all(sum(key)==n and tuple(sorted(key,reverse=True))==key for key in poly)
            assert helper.corner_value(poly)==F(corner['moments'][str(k)][n])
    for n,poly in enumerate(families['b']):assert all(sum(key)==n+1 for key in poly)
    portable_path=ROOT/'space_moments.json'
    if input_path!=portable_path.resolve():
        portable={'families':{name:[{','.join(map(str,key)):str(value) for key,value in sorted(poly.items())} for poly in families[name]] for name in names},
                  'raw_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(),
                  'backend_source_sha256':data['backend_source_sha256'],
                  'canonical_moment_count':data['canonical_moment_count'],
                  'production_seconds':data['seconds']}
        portable_path.write_text(json.dumps(portable,indent=2)+'\n')
    eta=F(corner['eta']);a=F(corner['a_selected']);cL=F(corner['c_lower'])
    assert eta==F(1,8) and a==F(679,1000)
    T=F(corner['T_lower']);coeff=list(map(F,corner['normalized_d_coefficients']))
    power={}
    for n in range(8):
        factor=F(11,4)*a**7*comb(7,n)*eta**n*T
        for k in range(1,5):
            for key,value in families['d'+str(k)][n].items():power[key]=power.get(key,F(0))+factor*coeff[k]*value
    for n in range(12):
        factor=-F(7,4)*a**11*comb(11,n)*eta**n
        for key,value in families['b'][n].items():power[key]=power.get(key,F(0))+factor*value
    power={k:v for k,v in power.items() if v}
    assert helper.corner_value(power)==F(corner['corner_polynomial_value'])
    assert max(sum(k) for k in power)==12 and max(max(k) for k in power)==6
    face,degrees,tensor,index=helper.bernstein_face(power)
    assert degrees==(6,6,6,6) and len(tensor)==2401
    minimum=tensor[index];lower=cL*minimum
    simple=F(lower.numerator*10**8//lower.denominator,10**8)
    assert all(value>=minimum for value in tensor.values())
    tensor_text=json.dumps({','.join(map(str,k)):str(v) for k,v in sorted(tensor.items())},sort_keys=True)
    output={'status':'one_registered_K3_full_space_candidate','new_C5_upper_bound':False,
            'positive_global_certificate':minimum>0,'eta':str(eta),'a':str(a),'c_lower':str(cL),
            'symmetric_power_coefficients':{','.join(map(str,k)):str(v) for k,v in sorted(power.items())},
            'partition_coefficient_count':len(power),'face_power_coefficient_count':len(face),
            'face_Bernstein_degrees':degrees,'Bernstein_coefficient_count':len(tensor),
            'minimum_Bernstein_index':index,'minimum_Bernstein_value':str(minimum),
            'global_g_exact_rational':str(lower) if minimum>0 else None,
            'global_g_simple_rational':str(simple) if minimum>0 else None,
            'strictly_improves_4583_over_1e8':minimum>0 and simple>F(4583,10**8),
            'float_display_only':float(lower),'Bernstein_coefficients_sha256':hashlib.sha256(tensor_text.encode()).hexdigest(),
            'crosschecks':{'all_004_space_coefficients_equal':True,'all_006_corner_moments_equal':True,
                           'new_d3_d4_corner_check_uses_same_graph_backend':True},
            'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'registration.md',ROOT/'graph_moments.cpp',Path(__file__),portable_path,corner_path,old_path,helper_path]},
            'runtime':{'python':sys.executable,'seconds':time.monotonic()-began}}
    (ROOT/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ['symmetric_power_coefficients','source_sha256']},indent=2))


if __name__=='__main__':main()
