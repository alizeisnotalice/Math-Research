#!/usr/bin/env python3
"""Odd-dimensional ball translation loss: exact polynomial slice integration."""
from fractions import Fraction as F
from math import comb,factorial
from pathlib import Path
import argparse,json

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    batches=[]
    for ns in [(1,3,5),(9,17),(33,65)]:
        rows=[]
        for n in ns:
            m=(n-1)//2;beta=F(factorial(2*m+1),2**(2*m+1)*factorial(m)**2)
            for q in [F(1,16),F(1,4),F(1,2),F(1)]:
                loss=2*beta*sum(((-1)**i*comb(m,i)*(q/2)**(2*i+1)/F(2*i+1) for i in range(m+1)),F(0))
                assert 0<=loss<=min(F(1),beta*q)
                rows.append(dict(n=n,shift_over_radius=str(q),beta=str(beta),one_sided_loss=str(loss),bound=str(min(F(1),beta*q))))
        batches.append(rows)
    args.output.write_text(json.dumps(dict(status='passed',arithmetic='Fraction polynomial integral',batches=batches,case_count=sum(map(len,batches))),indent=2)+'\n')
    print('28 exact slice checks passed')
if __name__=='__main__':main()
