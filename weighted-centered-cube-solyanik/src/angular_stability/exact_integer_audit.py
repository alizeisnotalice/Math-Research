#!/usr/bin/env python3
"""Independent Fraction audit of scalar, averaged, and sharp-sequence rows."""
from fractions import Fraction as Q
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
PS=(2,3,4,8,16)

def balanced(d): return [Q(1,2),Q(1,2)],[1-d,1+d]
def hole(q): return [q,1-q],[q,1+q]
def spike(p,scale):
    M=Q(2**scale); eps=M**(-(p+1)); c=(1-eps*M)/(1-eps)
    return [eps,1-eps],[M,c]
def random_vec(seed):
    w0=[Q((seed*7+3*i)%11+1,37) for i in range(4)]; w=[z/sum(w0,Q(0)) for z in w0]
    x0=[Q((seed*13+5*i)%29+1,(seed+i)%7+2) for i in range(4)]
    mean=sum((a*b for a,b in zip(w,x0)),Q(0)); return w,[z/mean for z in x0]

def check(w,x,p):
    assert sum((a*b for a,b in zip(w,x)),Q(0))==1
    B=lambda z:z**p-1-p*(z-1)
    scalar=all(abs(z-1)**p<=B(z) and (p-1)*(z-1)**2<=B(z) for z in x)
    b=sum((a*B(z) for a,z in zip(w,x)),Q(0)); dp=sum((a*abs(z-1)**p for a,z in zip(w,x)),Q(0)); d2=sum((a*(z-1)**2 for a,z in zip(w,x)),Q(0))
    if p==2: valid=scalar and dp==b and d2==b
    else: valid=scalar and dp<=b and (p-1)*d2<=b
    return valid,b,dp,d2

def main():
    count=0; failed=[]; sharp=[]
    for p in PS:
        family=[]
        family += [('balanced',*balanced(Q(1,2**j))) for j in (3,8,16)]
        family += [('hole',*hole(Q(1,2**j))) for j in (3,6,9)]
        family += [(f'spike_{j}',*spike(p,j)) for j in (1,3,10)]
        family += [('random',*random_vec(s)) for s in range(1,16)]
        for label,w,x in family:
            valid,b,dp,d2=check(w,x,p);count+=1
            if not valid: failed.append({'p':p,'label':label,'weights':[str(z) for z in w],'values':[str(z) for z in x]})
            if label.startswith('spike'): sharp.append({'p':p,'label':label,'ratio_Dp_over_B':str(dp/b)})
            if label=='hole': sharp.append({'p':p,'label':'hole','ratio_(p-1)D2_over_B':str((p-1)*d2/b)})
    out={'scope':'exact Fraction scalar and averaged checks for integer p only','cases':count,
         'failed':failed,'sharp_sequence_ratios':sharp}
    (ROOT/'exact_integer_audit.json').write_text(json.dumps(out,indent=2))
    assert not failed
    print('exact_integer_cases=%d failed=0 sharp_rows=%d' % (count,len(sharp)))

if __name__=='__main__': main()
