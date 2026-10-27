import csv,json,math,platform,sys,time
from fractions import Fraction as Q
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from verify import certify,ceilq
from pathlib import Path
p=Path(__file__).parent
rows=[]; trace=[]
for n in list(range(2,21))+[32,64]:
    calls=[]
    def f(y):
        t=math.exp(y); h=(1+1/t)/2
        integ=quad(lambda s: 1-0.5*(1+1/(2*t)-0.5*math.exp(-s))**(n-1),0,y,epsabs=1e-11)[0]
        b=1-h**n+n*integ; a=(1+1/t)**n
        value=(math.sqrt(a)+math.sqrt(b))**2
        calls.append([y,t,value]); return value
    if n==7:
        t=Q(3391749,524288); f(math.log(float(t))); source='original seven-dimensional certificate'
    else:
        opt=minimize_scalar(f,bounds=(1e-9,20),method='bounded',options={'maxiter':94,'xatol':1e-10})
        raw=math.exp(opt.x); t=Q(round(raw*2**20),2**20)
        if t<=1: t=Q(2**20+1,2**20)
        source='bounded double candidate rounded to denominator 2^20'
    assert len(calls)<=96
    A,lo,hi,fl,fu=certify(n,t)
    U=Q(4910992105193,200000000000) if n==7 else ceilq(fu,10**10)+Q(1,10**10)
    assert fu<U
    row={'n':n,'t':str(t),'source':source,'A':str(A),'B_lower':str(lo),'B_upper':str(hi),'F_lower':str(fl),'F_upper':str(fu),'U':str(U),'boundary_2_pow_n':str(2**n),'candidate_F_display':format(float((fl+fu)/2),'.12g'),'endpoint_comparison':'candidate below endpoint' if fu<2**n else 'endpoint is smaller or not separated','evaluations':len(calls)}
    rows.append(row); trace.append({'n':n,'evaluations':calls,'selected_t':str(t),'selection':source})
    print(n,str(t),row['candidate_F_display'],row['endpoint_comparison'],len(calls),flush=True)
(p/'table.json').write_text(json.dumps(rows,indent=2)+'\n')
(p/'trace.json').write_text(json.dumps(trace,indent=2)+'\n')
with (p/'table.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
(p/'table.tex').write_text('\n'.join(f"{r['n']} & $\\frac{{{Q(r['t']).numerator}}}{{{Q(r['t']).denominator}}}$ & {r['candidate_F_display']} & ${r['U']}$ & ${r['boundary_2_pow_n']}$ \\\\" for r in rows)+'\n')
(p/'environment.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'scipy':__import__('scipy').__version__,'method':'scipy bounded search in log(t), positive integral; rational exact verification','max_evaluations':96,'note':'Search is heuristic and does not establish global optimality. U is a strict upper bound.'},indent=2)+'\n')
# Add finite-a certificates after writing the selected points; no new search.
import add_finite
