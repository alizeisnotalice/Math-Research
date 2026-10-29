"""Independent exact geometry/metadata checks of new general-R source. No cache or MC."""
from pathlib import Path
from fractions import Fraction as F
import sys,importlib.util,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;OTHER=HERE.parent/'round040_source'
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('new40',OTHER/'source_tables.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
validation=json.loads((OTHER/'SOURCE_VALIDATION.json').read_text())
assert validation['source_SHA256']==hashlib.sha256((OTHER/'source_tables.py').read_bytes()).hexdigest()
assert validation['validation_SHA256']==hashlib.sha256((OTHER/'validate_source.py').read_bytes()).hexdigest()
cases=[]
for R in [4,8,16,32]:
    src=mod.Source(R);half=F(1,2*src.w*src.v);touch=0;tests=0;mingap=F(1)
    for b in range(R):
        positions=[]
        for j,label in enumerate(src.labels):
            expected=(F(int(label),4)+F(2*b+1 if j<src.cut else 1-2*b,2*R))%1
            actual=F(int(src.atom_offset_numerators[b,j]),2*R)
            assert expected==actual
            assert 0<=actual-half<=actual+half<=1
            touch+=(actual-half==0)+(actual+half==1);tests+=1
            positions.append((F(j)+actual)/src.K)
        for group in [positions[:src.cut],positions[src.cut:]]:
            mingap=min(mingap,*[z-a for a,z in zip(group,group[1:])],group[0]+1-group[-1])
    assert mingap>=F(1,4*src.K)>F(1,src.w*src.K*src.v)+F(1,src.Q)
    assert src.background_lower==float(F(437,918))
    cases.append({'R':R,'exact_offset_and_box_tests':tests,'touching_edge_count':touch,
                  'minimum_within_class_circular_gap':str(mingap),'source_mass_exact':'1',
                  'common_plus_mass_formula':'127/255 + (1/1020) sum c_m sin(2pi m 127/255)',
                  'positive_background_rational':str(F(437,918))})
out={'cases':cases,'frozen_source_hash_matches_validation':True,'validation_hash_matches':True,
     'scope':'Exact geometry of ideal rational boxes; D/eta actual FP remains un-certified; no large cache/MC'}
p=HERE/'source_checks.json'
if p.exists():raise RuntimeError('Refusing replacement')
p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
