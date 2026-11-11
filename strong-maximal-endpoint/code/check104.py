"""Exact finite envelope computation for Boolean subcubes."""
import argparse,math,json
from fractions import Fraction
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def bounds(d):
    mass=2**d
    coeff=1
    upper=Fraction(0);lower=Fraction(0);excluded=0
    for s in range(d+1):
        a=Fraction(coeff,2**s)
        b=Fraction((2**(d-s))*coeff,mass-excluded)
        lower=max(lower,a);upper=max(upper,b)
        excluded+=coeff
        coeff=coeff*(d-s)//(s+1)
    assert lower<=upper
    return {'dimension':d,'lower_num':str(lower.numerator),'lower_den':str(lower.denominator),'upper_num':str(upper.numerator),'upper_den':str(upper.denominator),'relative_width':float(upper/lower-1)}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('output',type=Path);a=ap.parse_args();ds=list(range(1,257))+list(range(320,2049,64))
    with ProcessPoolExecutor(12) as pool:rows=list(pool.map(bounds,ds))
    a.output.write_text(json.dumps({'cases':len(rows),'selected':[r for r in rows if r['dimension'] in (2,3,4,5,16,64,256,1024,2048)],'max_dimension':2048},indent=2)+'\n')
