from pathlib import Path
import importlib.util,json,time,signal
root=Path(__file__).resolve().parent;out=root/'union_cases_reproduced.json';assert not out.exists()
def stop(*args):raise TimeoutError('10s fixed six-case budget')
signal.signal(signal.SIGALRM,stop);signal.alarm(10);t=time.monotonic()
up=root
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
one=load('first_union',up/'check_exact2.py');two=load('second_union',up/'lb_check2.py')
cases=[[],[(0,1,0,1),(1,2,0,1)],[(0,2,0,2),(1,3,1,3)],[(0,3,0,3),(1,2,1,2)],[(0,2,0,2)]*3,[(0,0,0,2),(0,2,1,1),(0,1,0,1)]]
rows=[]
for boxes in cases:
 exact=sum(any(a<=x<b and c<=y<d for a,b,c,d in boxes) for x in range(3) for y in range(3))
 a,b=one.union_area_int(boxes),two.union_area(boxes);assert a==b==exact
 rows.append({'rectangles':boxes,'integer_grid_area':exact,'first':a,'second':b})
signal.alarm(0);out.write_text(json.dumps({'status':'passed','cases':rows,'seconds':time.monotonic()-t,'budget_seconds':10,'scope':'Six fixed union cases only; not a configuration certificate.'},indent=2)+'\n');print(out.read_text())
