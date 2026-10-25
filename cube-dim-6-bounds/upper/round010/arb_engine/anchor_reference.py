"""Serialize exact outward endpoints of the same six independent A007 DD anchors."""
from pathlib import Path
from fractions import Fraction
import importlib.util,json,hashlib
from flint import arb,ctx
HERE=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper/round007/g_time_moments/run_moments.py')
spec=importlib.util.spec_from_file_location('A007_reference',SRC)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);ctx.prec=640
reg=json.loads((HERE/'registration.json').read_text())
native=json.loads((HERE/'anchor_native.json').read_text())

def exact_bounds(row):
    e=int(row['common_binary_exponent'])
    return [Fraction(int(row[k]))*Fraction(2)**e for k in ['lower_mantissa','upper_mantissa']]

rows=[]
for i,(nodes,n)in enumerate(zip(reg['registered_anchors'],native['anchors'])):
    r=mod.dd([tuple(v)for v in nodes],arb(1)/250)
    lm,le=map(int,r.lower().man_exp());um,ue=map(int,r.upper().man_exp());e=min(le,ue)
    endpoint={'lower_mantissa':str(lm<<(le-e)),'upper_mantissa':str(um<<(ue-e)),
              'common_binary_exponent':str(e)}
    lo,hi=exact_bounds(n);rl,ru=exact_bounds(endpoint)
    assert lo<=rl<=ru<=hi
    rows.append({'index':i,'native_contains_python640_reference_interval':True,
                 'python640_exact_outward_dyadic_interval':endpoint})
result={'status':'PASS_ALL6','reference_source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),
        'precision_python_ref':640,'precision_native':512,'native_version':native['flint_version'],
        'native_sizeof_arb':native['sizeof_arb'],'same_six_registered_anchor_cases':True,'anchors':rows}
(HERE/'anchor_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS six same anchors: exact dyadic reference endpoints retained')
