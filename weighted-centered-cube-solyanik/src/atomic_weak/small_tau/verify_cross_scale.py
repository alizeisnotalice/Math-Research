#!/usr/bin/env python3
"""Independent all-subsets audit for each best corrected small-tau candidate."""
import json
import argparse
from pathlib import Path
from decimal import Decimal as D

import small_tau_atomic_search as core


def merge(intervals):
    answer=[]
    for left,right in sorted(intervals,key=lambda row:row[0]):
        if left>=right: continue
        if answer and left<=answer[-1][1]: answer[-1]=(answer[-1][0],max(answer[-1][1],right))
        else: answer.append((left,right))
    return answer


def all_subsets_halo(atoms):
    intervals=[]; n=len(atoms)
    for mask in range(1,1<<n):
        selected=[atoms[j] for j in range(n) if mask>>j&1]
        mass=sum((m for x,m in selected),D(0))
        intervals.append((max(x for x,m in selected)-mass/2,
                          min(x for x,m in selected)+mass/2))
    return merge(intervals)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,default=Path(__file__).resolve().parent/'final_candidates.json')
    records=json.loads(parser.parse_args().input.read_text())
    for tau_text in ('1/64','1/128','1/256'):
        row=max((r for r in records if r['tau']==tau_text),key=lambda r:D(r['tau_times_score']))
        atoms=[(D(x),D(m)) for x,m in row['atoms']]
        assert core.block_halo(atoms)==all_subsets_halo(atoms)
        tau=core.parse_tau(tau_text)
        value=tau*core.score(atoms,tau)
        assert abs(value-D(row['tau_times_score'])) <= D('1e-70')
        print(tau_text,row['family'],row['atom_count'],value,'ALL-SUBSETS-PASS')


if __name__=='__main__': main()
