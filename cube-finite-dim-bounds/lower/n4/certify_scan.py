"""Exact rational outer covers of this batch's finite witness unions, not family bounds."""
from fractions import Fraction as F
from itertools import product
from functools import lru_cache
from pathlib import Path
import json, hashlib, time
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[1]
H=10**12

def add(x,y):return (x[0]+y[0],x[1]+y[1])
def scale(x,t):return (x[0]*t,x[1]*t) if t>=0 else (x[1]*t,x[0]*t)
def mul(x,y):
    z=[a*b for a,b in product(x,y)];return min(z),max(z)
def floor(x):return x.numerator//x.denominator
def ceil(x):return -floor(-x)
def rootfloor(x):
    lo=0;hi=1
    while hi**4<=x:hi*=2
    while lo+1<hi:
        mid=(lo+hi)//2
        if mid**4<=x:lo=mid
        else:hi=mid
    return lo
cert=json.loads((ROOT/'baseline/original_certificates.json').read_text())[1]
D=cert['root_denominator'];nums=cert['root_lower_numerators']
for num,target in zip(nums,[2,3]):assert num**4<target*D**4<(num+1)**4
u=(F(nums[0],D),F(nums[0]+1,D));v=(F(nums[1],D),F(nums[1]+1,D))
A=[add((F(1),F(1)),scale(mul(u,u),F(1,2))),scale(add((F(2),F(2)),mul(mul(mul(u,u),u),v)),F(1,4))]+[scale(add((F(1),F(1)),u),F(1,2))]*2

@lru_cache(None)
def radius(mass):
    assert mass>0
    dr=10**16;k=rootfloor(floor(mass*dr**4));return F(k,2*dr),F(k+1,2*dr)

def boxes(f,w,pattern):
    assert F(0)<=f<=F(1,2) and F(0)<=w<=F(7,10)
    assert pattern in ('single','diagonal')  # q components nonnegative; sm/sM use this.
    q=[scale(a,f*int(pattern=='diagonal' or j==0)) for j,a in enumerate(A[1:])]
    L=[scale(A[0],2)]+A[1:];rows=set()
    for ell in product(range(1,5),repeat=4):
        for k1 in range(-4,5):
            parity=1 if k1%2==0 else -1
            mass=F(ell[1]*ell[2]*ell[3])*(ell[0]+w*(parity if ell[0]%2 else 0))
            r=radius(mass)
            lo1=add(scale(A[0],k1+ell[0]-1),scale(r,-1))[0]
            hi1=add(scale(A[0],k1),r)[1]
            if hi1<=0 or lo1>=L[0][1] or hi1<=lo1:continue
            sm=[scale(x,F(parity,2)) for x in q] if ell[0]==1 else [scale(x,F(-1,2)) for x in q]
            sM=[scale(x,F(parity,2)) for x in q] if ell[0]==1 else [scale(x,F(1,2)) for x in q]
            choices=[]
            for j in range(1,4):
                opts=[]
                for k in range(-4,5):
                    lo=add(add(scale(A[j],k+ell[j]-1),sM[j-1]),scale(r,-1))[0]
                    hi=add(add(scale(A[j],k),sm[j-1]),r)[1]
                    lo=max(F(0),lo);hi=min(L[j][1],hi)
                    if hi>lo:opts.append((floor(lo*H),ceil(hi*H)))
                choices.append(opts)
            for zz in product(*choices):rows.add((floor(max(F(0),lo1)*H),ceil(min(L[0][1],hi1)*H))+sum(zz,()))
    b=np.array(sorted(rows),dtype=np.int64);lo=b[:,::2];hi=b[:,1::2];keep=[]
    for i in range(len(b)):
        contained=np.all(lo<=lo[i],axis=1)&np.all(hi>=hi[i],axis=1);contained[i]=False
        if not contained.any():keep.append(i)
    return tuple(tuple(map(int,b[i])) for i in keep)

def volume(box):
    @lru_cache(None)
    def rec(B):
        if not B:return 0
        if len(B[0])==2:
            a,b=B[0];s=0
            for c,e in B[1:]:
                if c>b:s+=b-a;a,b=c,e
                else:b=max(b,e)
            return s+b-a
        cuts=sorted({x for b in B for x in b[:2]})
        return sum((v-u)*rec(tuple(sorted({b[2:] for b in B if b[0]<=u and b[1]>=v}))) for u,v in zip(cuts,cuts[1:]))
    return rec(tuple(sorted(box)))

def subset_boxes(row,inner=False):
    f=F(str(row['f']));w=F(str(row['w']));pattern=row['pattern']
    assert F(0)<=f<=F(1,2) and F(0)<=w<=F(7,10) and pattern in ('single','diagonal')
    q=[scale(a,f*int(pattern=='diagonal' or j==0)) for j,a in enumerate(A[1:])]
    L=[scale(A[0],2)]+A[1:]
    ids=list(product(range(-3,5),repeat=4));out=[]
    for subset in row['subsets']:
        coords=[];mass=F(0)
        for ix in subset:
            k=ids[ix];parity=1 if k[0]%2==0 else -1
            mass+=1+w*parity
            coords.append([scale(A[0],k[0])]+[add(scale(A[j],k[j]),scale(q[j-1],F(parity,2))) for j in range(1,4)])
        r=radius(mass);box=[]
        for j in range(4):
            lo=max(F(0),max(p[j][1 if inner else 0] for p in coords)-r[0 if inner else 1])
            hi=min(L[j][0 if inner else 1],min(p[j][0 if inner else 1] for p in coords)+r[0 if inner else 1])
            box.extend([ceil(lo*H),floor(hi*H)] if inner else [floor(lo*H),ceil(hi*H)])
        if all(box[2*j+1]>box[2*j] for j in range(4)):out.append(tuple(box))
    return tuple(out)

def check_index_margin():
    """Exact labels, independent of geometric volume; paired labels need not lie inside F."""
    margin=8
    def label(k):return (k[0]//2,k[1],k[2],k[3])
    per_axis={start+offset for start in range(-4,5) for ell in range(1,5) for offset in range(ell)}
    assert per_axis==set(range(-4,8))
    rectangular_atoms=list(product(sorted(per_axis),repeat=4))
    rectangular_labels=[label(k) for k in rectangular_atoms]
    assert all(max(map(abs,h))<=margin for h in rectangular_labels)
    ids=list(product(range(-3,5),repeat=4))
    rows=json.loads((P/'captured_subset_results.json').read_text())['rows']
    extra_indices=[ids[ix] for row in rows for subset in row['subsets'] for ix in subset]
    extra_labels=[label(k) for k in extra_indices]
    assert all(max(map(abs,h))<=margin for h in extra_labels)
    # Exact translation identity for all represented first-coordinate indices and cell labels.
    for k1 in range(-4,8):
        for cell1 in range(-9,10):assert (k1+2*cell1)//2==k1//2+cell1
    for w in map(lambda x:F(str(x)),[0,.02,.08,.2,.4,.7]):assert (1+w)+(1-w)==2
    def extrema(labels):return [[min(h[j] for h in labels),max(h[j] for h in labels)] for j in range(4)]
    return dict(margin=margin,label_map=['floor(k1/2)','k2','k3','k4'],
        paired_representatives=['(2h1,h2,h3,h4)','(2h1+1,h2,h3,h4)'],
        pair_mass='(1+w)+(1-w)=2',rectangular_original_index_range=[-4,7],
        rectangular_cover_atom_count=len(rectangular_atoms),rectangular_label_extrema=extrema(rectangular_labels),
        added_subset_count=sum(len(row['subsets']) for row in rows),added_atom_occurrences=len(extra_indices),
        added_subset_label_extrema=extrema(extra_labels),reception_cell_labels='[-N,N]^4',
        retained_paired_labels='[-N-8,N+8]^4',retained_mass='2*(2N+17)^4',
        finite_weak_ratio_lower='(A_tau/2)*((2N+1)/(2N+17))^4',
        translation_identity='label(k+(2s1,s2,s3,s4))=label(k)+s',
        all_index_checks_passed=True)

def main():
    st=time.monotonic();data=json.loads((P/'staggered_results.json').read_text());oldlo=F(cert['lower_bound_interval'][0]);out=[]
    for row in data['rows']:
        if row['f']==row['w']==0:continue
        b=boxes(F(str(row['f'])),F(str(row['w'])),row['pattern']);up=F(volume(b),2*H**4)
        assert up<oldlo
        out.append(dict(f=str(row['f']),w=str(row['w']),pattern=row['pattern'],outer_ratio_upper=str(up),strictly_less_than_exact_old_lower=True,box_count=len(b)))
    result=dict(scope='Only finite witness unions: lengths 1..4 at registered parameters. No exclusion of other subsets or radii.',rounding_denominator=H,source_sha256=hashlib.sha256((ROOT/'baseline/original_certificates.json').read_bytes()).hexdigest(),old_exact_lower=str(oldlo),maximum_outer_upper=str(max(F(x['outer_ratio_upper']) for x in out)),elapsed=time.monotonic()-st,certified_cases=len(out),rows=out)
    enriched=[]
    for row in json.loads((P/'captured_subset_results.json').read_text())['rows']:
        base=boxes(F(str(row['f'])),F(str(row['w'])),row['pattern']);extra=subset_boxes(row)
        up=F(volume(base+extra),2*H**4);assert up<oldlo
        inner=subset_boxes(row,inner=True)
        gain=F(volume(base+inner)-volume(base),2*H**4)
        if extra:assert gain>0
        enriched.append(dict(strict_finite_enrichment_gain_lower=str(gain),f=str(row['f']),w=str(row['w']),pattern=row['pattern'],extra_subset_count=len(extra),outer_ratio_upper=str(up),strictly_less_than_exact_old_lower=True))
    result['enriched_rows']=enriched
    result['finite_block_transfer']=check_index_margin()
    result['input_sha256']={name:hashlib.sha256((P/name).read_bytes()).hexdigest() for name in ['preregister.md','staggered_results.json','captured_subset_results.json','nonrectangular_example.json']}
    result['elapsed']=time.monotonic()-st
    (P/'scan_certificate.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
if __name__=='__main__':main()
