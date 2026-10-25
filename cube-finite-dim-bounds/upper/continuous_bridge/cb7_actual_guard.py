"""Exact outer-sphere guard for two author-supplied fixed decimal inputs.
Shares the finite candidate theorem; not independent mathematical audit.
"""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
OUT=Path(__file__).resolve().parent
ROWS=[(3,'0.4089479272201547','0.005004477192161995','0.0009386738283819751',Q(16,5)),
      (4,'0.18157763153667977','0.004','0.0005',Q(18,5))]

def run():
    results=[]
    for n,Atext,B4text,B6text,R in ROWS:
        A,B4,B6=map(Q,(Atext,B4text,B6text))
        d=B6*R**6;e=B4*R**4-Q(15,n+8)*d
        const=-Q(3,n+2)*B4*R**4+Q(30,(n+8)*(n+4))*d
        candidates=[]
        for support in range(1,n+1):
            candidates.append({'kind':'uniform','support':support,'value':d/Q(support**2)+e/Q(support)+const})
            if d==0:continue
            T=-2*e/(3*d)
            for p in range(1,support):
                q=support-p
                if p==q:continue # any stationary family already has uniform value
                u=(1-q*T)/(p-q);v=T-u
                if u>0 and v>0:
                    candidates.append({'kind':'two_values','support':support,'p':p,'q':q,'u':str(u),'v':str(v),'value':d*(p*u**3+q*v**3)+e*(p*u**2+q*v**2)+const})
        selected=min(candidates,key=lambda x:x['value'])
        gmin=selected['value'];margin=A-R**(2-n)+gmin
        assert margin>0 and R**2>n
        for candidate in candidates:candidate['value']=str(candidate['value'])
        results.append({'n':n,'A_exact':str(A),'B4_exact':str(B4),'B6_exact':str(B6),'R_exact':str(R),
                        'min_G_exact':str(gmin),'max_normalized_F_exact':str(-margin),'strict_negative_margin_exact':str(margin),
                        'strict_negative_margin_float_display':float(margin),'selected_pattern':selected,'all_candidates':candidates,
                        'scope':'exact full-sphere negativity only; cube minimum and mass/cost not certified here'})
    source=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-upper/cube-finite-dim-bounds/upper/u06_results.json')
    rec={'status':'PASS two exact outer guards','input_rule':'supplied decimal strings interpreted as exact rational values; no parameters changed','results':results,
         'author_results_read_path':str(source),'author_results_read_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'independence':'exact implementation shares analytic enumeration/formulas; requires independent audit'}
    (OUT/'cb7_actual_guard_results.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
if __name__=='__main__':run()
