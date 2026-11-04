#!/usr/bin/env python3
"""Read-only audit of the three historical multiscale searches.

The printed comparison is floating-point exploratory evidence only.  Exact
Fraction controls verify that the contiguous-block halo representation equals
the all-subsets representation in small rational instances.
"""
import json
import math
import random
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parent
HISTORICAL=ROOT/'history'/'pre_signature_20261105'


def merge(intervals):
    answer=[]
    for a,b in sorted(intervals):
        if a>=b: continue
        if answer and a<=answer[-1][1]: answer[-1]=(answer[-1][0],max(answer[-1][1],b))
        else: answer.append((a,b))
    return answer


def blocks(atoms):
    atoms=sorted(atoms); output=[]
    for left in range(len(atoms)):
        total=F(0)
        for right in range(left,len(atoms)):
            total+=atoms[right][1]
            output.append((atoms[right][0]-total/2,atoms[left][0]+total/2))
    return merge(output)


def subsets(atoms):
    output=[]
    for mask in range(1,1<<len(atoms)):
        selected=[atoms[j] for j in range(len(atoms)) if mask>>j&1]
        total=sum((m for x,m in selected),F(0))
        output.append((max(x for x,m in selected)-total/2,
                       min(x for x,m in selected)+total/2))
    return merge(output)


def controls():
    rng=random.Random(20261105)
    for n in range(1,9):
        for _ in range(25):
            x=sorted(rng.sample(range(-40,41),n))
            atoms=[(F(v,17),F(rng.randint(1,20),19)) for v in x]
            assert blocks(atoms)==subsets(atoms)
    print('exact geometry controls: 200 rational all-subsets comparisons PASS')


def main():
    controls()
    baseline={0.25:3*(5/3)**0.25,0.5:1.5*math.sqrt(5/3)}
    records=json.loads((ROOT/'historical_candidates.json').read_text())['records']
    total=sum(z['evaluations'] for z in records)
    assert len(records)==32 and total==7083648
    print('family,tau,atoms,evaluations,best,baseline,gap')
    for row in records:
        n=len(row['atoms']); family=row['family']
        dim=(n//2+1) if family=='multiscale' else ((n+1) if family=='asymmetric' else 2*n-1)
        assert row['evaluations']==dim*row['popsize']*(row['maxiter']+2)
        assert 'Maximum number of iterations' in row['message']
        atoms=[(float(x),float(m)) for x,m in row['atoms']]
        halo=blocks(atoms)
        def primitive(x): return math.copysign(abs(x)**row['tau']/row['tau'],x)
        numerator=sum(primitive(b)-primitive(a) for a,b in halo)
        denominator=sum(m*abs(x)**(row['tau']-1) for x,m in atoms)
        value=numerator/denominator
        assert math.isclose(value,row['final_score'],rel_tol=1e-10,abs_tol=1e-10)
    for family in sorted({z['family'] for z in records}):
        for tau in (0.25,0.5):
            group=[z for z in records if z['family']==family and z['tau']==tau]
            best=max(group,key=lambda z:z['final_score'])
            assert best['final_score']<baseline[tau]
            print(f"{family},{tau},{len(best['atoms'])},{sum(z['evaluations'] for z in group)},"
                  f"{best['final_score']:.16f},{baseline[tau]:.16f},{baseline[tau]-best['final_score']:.3e}")
    print('historical new-search evaluations (both tau; excluding superseded free-v1 run):',total)


if __name__=='__main__': main()
