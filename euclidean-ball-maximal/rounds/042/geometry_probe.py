#!/usr/bin/env python3
"""Exact checks for round42 geometry no-go lemmas; not a row counterexample."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'runtime/work/general_input_20261003'))
from stage24_clock_replacement_probe import observer


def checked_clock(atoms, x, expected):
    g={j:sum((w for y,w in atoms if abs(x-y)<F(2)**(2-j)),F(0))/(2*F(2)**(2-j))
       for j in range(1,expected[1]+1)}
    J=next(j for j in g if g[j]>F(1,8))
    K=next(j for j in g if g[j]>F(1,2))
    assert (J,K)==expected
    return J,K


def h(cells,radii,y,k):
    return sum((max(F(0),min(F(v['hi']),y+radii[k])-max(F(v['lo']),y-radii[k]))/(2*radii[k])
                for v in cells if v['K']==k),F(0))


def main():
    shell=[]
    for n in (1,2,16,128):
        epsilon=F(1,1024*n)
        s=(1-epsilon/2)**n
        m=(1+s)/4
        mx=(1+s)/(2*s)
        assert s>F(1,3) and 1<mx<2
        assert m/4<F(1,8)<m/2 and m<F(1,2)<2*m
        shell.append(dict(n=n,epsilon=str(epsilon),near_mass=str(m),MP_at_origin=str(mx),
                          normalized_shell_mass=str(2*m/(n*epsilon)),J=2,K=4))
    lens=[]
    for k,tau in ((8,F(1,32)),(8,F(1,128)),(12,F(1,512))):
        r=F(2)**(2-k);R=F(1,4);y=R-tau*r;z=R+r-2*tau*r
        p=2*r*(1-tau/2)
        atoms=[(F(-1,8),F(3,8)),(y,p),(F(10),F(5,8)-p)]
        assert sum(w for _,w in atoms)==1
        checked_clock(atoms,F(0),(2,4));checked_clock(atoms,z,(2,k))
        assert 1<(1-tau/2)/(1-tau)<2
        assert F(3,2)<2 and (F(3,8)+p)/(2*y)<1
        assert (F(3,8)+p)/(2*(z+F(1,8)))<1
        lo=max(-R,z-r);hi=min(R,z+r)
        mass=sum((w for a,w in atoms if lo<a<hi),F(0))
        assert hi-lo==2*tau*r and mass==p
        saved,_,radii,_=observer(atoms,F(1),'band',4)
        cells=saved['observer_cells']
        # Full actual observer reconstruction verifies open neighboring cells too.
        for x,clock in ((F(0),(2,4)),(z,(2,k))):
            covering=[v for v in cells if F(v['lo'])<=x<=F(v['hi'])]
            assert covering and all((v['J'],v['K'])==clock for v in covering)
        ks={v['K'] for v in cells}
        row=sum((w*h(cells,radii,a,4)*sum((h(cells,radii,a,l) for l in ks if l>4),F(0))
                 for a,w in atoms),F(0))
        X4=sum((F(v['hi'])-F(v['lo']) for v in cells if v['K']==4),F(0))
        assert row<=F(3,2)*X4
        lens.append(dict(k=k,tau=str(tau),fine_source_mass=str(p),lens_length=str(hi-lo),
                         mass_over_lens_length=str(mass/(hi-lo)),coarse_clock=[2,4],fine_clock=[2,k],
                         actual_row4=str(row),X4=str(X4),row4_over_X4=str(row/X4)))
    result=dict(status='passed',arithmetic='Fraction',shell=shell,lens=lens,
                aggregate_row_disproved=False,thin_source_shell_and_lens_density_strengthenings_disproved=True)
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    output=parser.parse_args().output
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='passed',shell_dimensions=[v['n'] for v in shell],lens=lens,aggregate_row_disproved=False),indent=2))

if __name__=='__main__':main()
