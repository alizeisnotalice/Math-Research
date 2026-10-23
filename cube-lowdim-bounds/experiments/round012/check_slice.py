from pathlib import Path
import importlib.util,json,math
source=Path(__file__).resolve().parent.parent/'round008/search.py'
spec=importlib.util.spec_from_file_location('prior_slice',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def h4(t):
    return 8*(2*t*t+11/3+(1+t)**3*math.log(t/(1+t))+(0 if t==1 else (1-t)**3*math.log(t/abs(1-t))))
rows=[]
for t in [.17,.43,.79,1.,1.23]:
    actual=m.slice_generator_numeric(4,t);closed=h4(t)
    rows.append({'t':t,'closed':closed,'direct':actual,'abs_error':abs(actual-closed)})
assert max(r['abs_error'] for r in rows)<5e-7
Path(__file__).with_name('slice_check.json').write_text(json.dumps(rows,indent=2)+'\n')
print(rows)
