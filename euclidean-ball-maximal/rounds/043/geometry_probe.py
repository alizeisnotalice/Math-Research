#!/usr/bin/env python3
"""Exact finite checks of core cutoff and actual-clock local drift obstruction."""
from fractions import Fraction as F
from pathlib import Path
import argparse,json,sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'runtime/work/general_input_20261003'))
from stage24_clock_replacement_probe import observer

def core(n):
    lo,hi=F(0),F(1)
    target=F(1,(n+1)**2)
    for _ in range(48):
        mid=(lo+hi)/2
        if mid**n<=target:lo=mid
        else:hi=mid
    assert lo**n<=target<=hi**n
    def count(p):
        d=max(0,p.denominator.bit_length()-p.numerator.bit_length())
        while F(1,2**(d+1))>p:d+=1
        while d and F(1,2**d)<=p:d-=1
        return d
    low=count((1-lo)**n);high=count((1-hi)**n)
    assert low==high
    cost=2*low*target
    assert cost<=1
    return dict(n=n,rho_interval=[str(lo),str(hi)],N=low,normalized_cost=str(cost))

def local(bits):
    eps=F(1,2**bits);R=F(1,4);r=F(1,64);tau=F(1,32)
    y=R-tau*r;z=R+(1-2*tau)*r;p=2*r*(1-tau/2);m=F(1,4)-p-eps
    atoms=sorted([(F(-1,8),2*eps),(F(1,10),m),(y,p),(F(10),F(3,4)-eps)])
    assert all(w>0 for _,w in atoms) and sum(w for _,w in atoms)==1
    result,_,radii,_=observer(atoms,F(1),'band',4)
    cells=result['observer_cells'];width=tau*r/8
    for center,K in [(F(0),4),(z,8)]:
        length=sum((max(F(0),min(F(c['hi']),center+width)-max(F(c['lo']),center-width)) for c in cells if c['J']==2 and c['K']==K),F(0))
        assert length==2*width
    mass=lambda center,rad:sum((w for a,w in atoms if abs(a-center)<rad),F(0))
    gap=(mass(F(0),R)-mass(z,R))/(2*R)
    assert gap==4*eps and mass(z,r)==p
    assert abs(y)<R and z>R and abs(z-y)<r
    assert all(mass(y,radii[k])==mass(z,radii[k]) for k in range(1,9))
    def h4(center):
        return sum((max(F(0),min(F(c['hi']),center+R)-max(F(c['lo']),center-R))/(2*R) for c in cells if c['K']==4),F(0))
    drift=h4(y)-h4(z)
    assert drift==F(31,1024)
    positive_kernel=p/(4*R*r)
    return dict(epsilon=str(eps),coarse_gap=str(gap),positive_kernel=str(positive_kernel),
        ratio_to_gap=str(positive_kernel/gap),fixed_rectangle_area=str((2*width)**2),
        actual_neighbor_clocks_verified=True,first_eight_density_vectors_equal=True,h4_drift=str(drift))

def spectral(N):
    atoms=[(F(5*i),F(1,N)) for i in range(N)]
    result,_,radii,_=observer(atoms,F(1),'band',4)
    cells=result['observer_cells']
    A=[sum((max(F(0),min(F(c['hi']),y+radii[c['K']])-max(F(c['lo']),y-radii[c['K']]))/(2*radii[c['K']]) for c in cells),F(0)) for y,w in atoms]
    M=sum(w*a for (y,w),a in zip(atoms,A))
    transition=[[F(0) for _ in atoms] for _ in atoms]
    for c in cells:
        lo,hi=F(c['lo']),F(c['hi']);mid=(lo+hi)/2;r=radii[c['K']];v=2*r
        hits=[abs(y-mid)<r for y,w in atoms]
        g=sum((w for (y,w),hit in zip(atoms,hits) if hit),F(0))/v
        for i in range(N):
            for j in range(N):
                if hits[i] and hits[j]:transition[i][j]+=(hi-lo)*atoms[j][1]/(A[i]*v*v*g)
    assert all(transition[i][j]==int(i==j) for i in range(N) for j in range(N))
    assert A==[F(1,2)]*N and M==F(1,2)
    return dict(N=N,actual_kernel='identity',stationary_weights=[str(w) for _,w in atoms],
        source_A='1/2',M='1/2',spectral_gap='0',arithmetic='Fraction')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    batches=[(1,2,3),(4,8,16,32),(64,128,256,512)]
    data=dict(status='passed',arithmetic='Fraction; exact rational root brackets',
              core_batches=[[core(n) for n in ns] for ns in batches],
              local_batches=[[local(b) for b in bs] for bs in [(11,12),(16,24),(32,48)]],
              spectral_batches=[spectral(N) for N in (2,4,8)],
              scope='finite algebra/clock checks; general lemmas proved separately in report',main_bound_proved=False)
    args.output.write_text(json.dumps(data,indent=2)+'\n')
    print('passed: 11 core dimensions, 6 local parameters, 3 identity-kernel inputs')
if __name__=='__main__':main()
