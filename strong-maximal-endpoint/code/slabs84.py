"""Exact rational central Irwin--Hall slab probabilities; no floating cancellation."""
import json,math,hashlib,sys,time,multiprocessing as mp
from fractions import Fraction as Q
from pathlib import Path

def job(z):
    n,p=z;d=n**(p-1);c=n*d//2+1;total=0;binom=1
    for k in range(n//2+1):
        if k:binom=binom*(n-k+1)//k
        total+=(-1 if k%2 else 1)*binom*(c-k*d)**n
    prob=Q(2*total,math.factorial(n)*d**n)-1
    lower=Q(2,n**p)
    if not prob>lower:raise ArithmeticError('central slab lower bound failed')
    den=10**50;lo=Q((prob*den).__floor__(),den);hi=lo+Q(1,den)
    nb=abs(prob.numerator).to_bytes((abs(prob.numerator).bit_length()+7)//8,'big');db=prob.denominator.to_bytes((prob.denominator.bit_length()+7)//8,'big')
    return dict(N=n,p=p,probability_interval=[str(lo),str(hi)],proved_probability_lower=str(lower),
      exact_probability_numerator_bits=prob.numerator.bit_length(),exact_probability_denominator_bits=prob.denominator.bit_length(),
      exact_probability_SHA256=hashlib.sha256(len(nb).to_bytes(8,'big')+nb+db).hexdigest(),pid=__import__('os').getpid())

def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);start=time.monotonic();jobs=[(2**k,p) for k in range(1,11) for p in [2,3]]
    with mp.get_context('spawn').Pool(12) as pool:rows=pool.map(job,jobs)
    r=dict(pass_all=True,workers=12,actual_processes=len({r['pid'] for r in rows}),seconds=time.monotonic()-start,arithmetic='exact integers and rational probabilities; outward50decimal interval only for compact storage',rows=rows)
    (out/'slabs84.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(pass_all=True,cases=len(rows),actual_processes=r['actual_processes'],seconds=r['seconds'])))
if __name__=='__main__':main()
