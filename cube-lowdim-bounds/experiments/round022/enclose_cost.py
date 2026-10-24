"""Exact rational cost enclosure for a rational truncated-power/spline certificate.
Does not check positivity. No imports from search, Arb, or public checkers.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse,json,hashlib

def enclose(path,digits=60):
    data=json.loads(Path(path).read_text())
    alpha,R,h,c=(Q(data[k]) for k in ('alpha','R','h','c'))
    assert 1<alpha<2 and R>0 and h>0 and c>0
    total=Q(0);seen=set()
    for key,value in data['a'].items():
        i,j=map(int,key.split(','));assert 0<=j<=i<=data['N'] and (i,j) not in seen
        seen.add((i,j));total+=(1 if i==j==0 else 4 if j==0 or i==j else 8)*Q(value)
    exponent=2-alpha;n=exponent.denominator;D=10**digits
    target=R**exponent.numerator;z=target*D**n;N=z.numerator//z.denominator
    lo,hi=0,1<<((N.bit_length()+n-1)//n)
    while lo+1<hi:
        m=(lo+hi)//2
        if m**n<=N:lo=m
        else:hi=m
    assert Q(lo,D)**n<=target<Q(lo+1,D)**n
    factor=c*alpha/(2-alpha);offset=h*h*total/2
    lower=factor*Q(lo,D)+offset;upper=factor*Q(lo+1,D)+offset
    return {'certificate_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'alpha':str(alpha),'R':str(R),'root_degree':n,'root_power':exponent.numerator,
            'root_scale':str(D),'root_floor':str(lo),'sum_coefficients':str(total),
            'cost_factor':str(factor),'cost_offset':str(offset),
            'lower':str(lower),'upper':str(upper),'width':str(upper-lower)}

def main():
    p=argparse.ArgumentParser();p.add_argument('certificate');p.add_argument('--baseline',required=True)
    p.add_argument('--out',required=True);args=p.parse_args()
    out=Path(args.out)
    if out.exists():raise SystemExit('Output exists; inspect instead of overwrite.')
    a=enclose(args.certificate);b=json.loads(Path(args.baseline).read_text())
    baseline=Q(b['new_cost_lower']);upper=Q(a['upper'])
    a.update(baseline_lower=str(baseline),strict_gap_lower=str(baseline-upper),
             improves_baseline=upper<baseline,below_user_target=upper<Q(359123,100000),
             scope='Cost only: C1/C2/C4 and analytic comparison required separately.')
    out.write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2))
if __name__=='__main__':main()
