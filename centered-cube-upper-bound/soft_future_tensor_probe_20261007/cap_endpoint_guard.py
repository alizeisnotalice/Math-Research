"""Extract the actual nested cap oracle and test sigma=1 endpoint handling."""
from pathlib import Path
import ast,importlib.util,json,math
import numpy as np
D=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('model',D/'model.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
tree=ast.parse((D/'run.py').read_text());solve=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='solve');node=next(n for n in solve.body if isinstance(n,ast.FunctionDef) and n.name=='cap_trial')
code=compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),str(D/'run.py'),'exec')
records=[]
for value,expect in [(.1,'NUMERICALLY_ENCLOSED_FULL_FUTURE_NOT_FLOAT_CERTIFIED'),(2.,'UNKNOWN')]:
    scope=dict(m=m,np=np,n=1,rd=3,sigmas=np.array([0.,1.]),LF=np.array([2.]),h=np.ones((3,2,1)),blo=np.full((3,2,1),value),bhi=np.full((3,2,1),value))
    exec(code,scope);res=scope['cap_trial'](1.)
    assert res['status']==expect and res['degenerate_sigma1_checked'] and res['tested_cells']==1,res
    records.append(dict(coordinate_soft_value=value,expected=expect,result=res))
(D/'cap_endpoint_validation.json').write_text(json.dumps(dict(status='PASS_ACTUAL_ORACLE_DEGENERATE_SIGMA1_LOGIC',records=records,scope='Synthetic coefficient check only; not a full-source FIRST or floating interval certificate.'),indent=2))
print('SIGMA1 CAP BRANCH PASS')
