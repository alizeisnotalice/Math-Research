"""Exact beta-mass enclosure and exterior-structure check; no positivity of the full generator.
Does not check positivity. No imports from search, Arb, or public checkers.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse,json,hashlib

def enclose(path,digits=60):
    data=json.loads(Path(path).read_text())
    alpha=Q(data['alpha']);beta=Q(data['beta'])
    R,h,c=(Q(data[k]) for k in ('R','h','c'))
    assert 1<alpha<=beta<2 and R>0 and h>0 and c>0
    total=Q(0);seen=set()
    for key,value in data['a'].items():
        i,j=map(int,key.split(','));assert 0<=j<=i<=data['N'] and (i,j) not in seen
        seen.add((i,j));total+=(1 if i==j==0 else 4 if j==0 or i==j else 8)*Q(value)
    exponent=2-beta;n=exponent.denominator;D=10**digits
    target=R**exponent.numerator;z=target*D**n;N=z.numerator//z.denominator
    lo,hi=0,1<<((N.bit_length()+n-1)//n)
    while lo+1<hi:
        m=(lo+hi)//2
        if m**n<=N:lo=m
        else:hi=m
    assert Q(lo,D)**n<=target<Q(lo+1,D)**n
    factor=c*beta/(2-beta);offset=h*h*total/2
    lower=factor*Q(lo,D)+offset;upper=factor*Q(lo+1,D)+offset
    return {'certificate_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'alpha':str(alpha),'beta':str(beta),'R':str(R),'root_degree':n,'root_power':exponent.numerator,
            'root_scale':str(D),'root_floor':str(lo),'sum_coefficients':str(total),
            'cost_factor':str(factor),'cost_offset':str(offset),
            'lower':str(lower),'upper':str(upper),'width':str(upper-lower)}

def structure(candidate_path,baseline_path):
    raw=Path(baseline_path).read_bytes();base=json.loads(raw)
    expected='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
    assert hashlib.sha256(raw).hexdigest()==expected
    d=json.loads(Path(candidate_path).read_text());w=d['structure']
    assert w['baseline_sha256']==expected
    assert (Q(d['alpha']),Q(d['R']),Q(d['h']),d['N'])==(Q(base['alpha']),Q(base['R']),Q(base['h']),base['N'])
    assert Q(d['c'])>0 and Q(w['s'])>=0 and Q(d['alpha'])<=Q(d['beta'])<2
    def parse(mapping):
        out={}
        for key,value in mapping.items():
            i,j=map(int,key.split(','));assert 0<=j<=i<=d['N'] and (i,j) not in out
            out[i,j]=Q(value)
        return out
    a,b,T=parse(d['a']),parse(base['a']),parse(w['T']);scale=Q(w['s'])
    assert all(i+j<=48 for (i,j),v in T.items() if v)
    assert all(a.get(ij,Q(0))==scale*b.get(ij,Q(0))+T.get(ij,Q(0)) for ij in set(a)|set(b)|set(T))
    support=Q(52)*Q(d['h']);assert support==Q(13,8)<Q(base['R'])
    return {'passed':True,'baseline_sha256':expected,'s':str(scale),'inner_support_bound':str(support),
            'T_nonzero_classes':sum(v!=0 for v in T.values()),'inherited_nonnegative_region':'r >= R',
            'transition_region':'13/8 <= r < R still requires full-kernel positivity; no c*mu >= s*c1 assertion is made',
            'scope':'K nonnegativity for r >= R inherits the already certified K1; all r < R positivity/dominance and global generator still require certification.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('certificate');p.add_argument('--baseline',required=True)
    p.add_argument('--kernel-baseline',required=True);p.add_argument('--out',required=True);args=p.parse_args()
    out=Path(args.out)
    if out.exists():raise SystemExit('Output exists; inspect instead of overwrite.')
    struct=structure(args.certificate,args.kernel_baseline)
    a=enclose(args.certificate);a['structure']=struct
    baseline_raw=Path(args.baseline).read_bytes()
    assert hashlib.sha256(baseline_raw).hexdigest()=='1f7b1a0184101a93f5d64e463377ee4e86374660fee3e3c2de7676ff46237159'
    a['baseline_arithmetic_sha256']=hashlib.sha256(baseline_raw).hexdigest();b=json.loads(baseline_raw)
    baseline=Q(b['new_cost_lower']);upper=Q(a['upper'])
    a.update(baseline_lower=str(baseline),strict_gap_lower=str(baseline-upper),
             improves_baseline=upper<baseline,below_user_target=upper<Q(17,5), user_target="17/5",
             scope='Cost only: C1/C2/C4 and analytic comparison required separately.')
    out.write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2))
if __name__=='__main__':main()
