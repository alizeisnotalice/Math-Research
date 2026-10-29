#!/usr/bin/env python3
"""Independent exact trace matrix and label-map rank (no q-state enumeration)."""
import pathlib,json
from certify import mul,forms
HERE=pathlib.Path(__file__).resolve().parent

def trace(x,d,p):
    z=0
    for _ in range(d):z^=x;x=mul(x,x,p)
    assert z in [0,1]
    return z

def rank(rows):
    piv={}
    for x in rows:
        while x:
            j=x.bit_length()-1
            if j in piv:x^=piv[j]
            else:piv[j]=x;break
    return len(piv)

def main():
    m=json.loads((HERE/'MANIFEST.json').read_text());d=m['d'];p=m['modulus'];g=m['generator'];basis=m['basis']
    def row(x):return sum(trace(mul(x,1<<j,p),d,p)<<j for j in range(d))
    pair=[row(1<<i) for i in range(d)];beta=[row(x) for x in basis]
    W=forms(basis);examples=[]
    # Generic lags 1 and 2: forbidden short-ratio exclusion gives zero intersection.
    for k,t in [(1,g),(2,mul(g,g,p))]:
        r=rank(basis+[mul(t,x,p) for x in basis]);assert r==2*m['b']
        examples.append({'k':k,'scalar':t,'joint_trace_rank':r,'expected':2*m['b']})
    # A concrete defect ratio W[1]/W[2], calculated by independent polynomial arithmetic.
    from certify import power
    t=mul(W[1],power(W[2],m['K']-1,p),p)
    r=rank(basis+[mul(t,x,p) for x in basis]);assert r==2*m['b']-1
    assert rank(pair)==d and rank(beta)==m['b']
    out={'trace_pairing_rows':pair,'trace_pairing_rank':rank(pair),'label_rows':beta,'label_map_rank':rank(beta),'short_lag_joint_ranks':examples,'sample_defect':{'a':1,'c':2,'scalar_ratio':t,'joint_trace_rank':r},'exact_full_field_label_fiber_count':m['q']//m['B'],'exact_nonzero_field_zero_label_fiber_count':m['q']//m['B']-1,'scope':'Linear algebra certificate; no full-period label list.'}
    (HERE/'TRACE_CHECK.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
