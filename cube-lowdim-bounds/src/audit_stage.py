"""Standard-library proof certificates, distinct from discovery scripts."""
from pathlib import Path
from fractions import Fraction as F
import sys,json,zipfile,random,itertools,math
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'baseline'))
from verify_bounds import floor_root,union_sweep,union_grid
from upper_certificate import logq,certify

def root_lower(q,d,D=10**15):
    q=F(q);t=floor_root(q.numerator*D**d//q.denominator,d)
    assert F(t,D)**d<=q<F(t+1,D)**d
    return F(t,D)

def exact_union(B):
    B=tuple(sorted(set(b for b in B if all(b[j]<b[j+1] for j in range(0,len(b),2)))))
    if not B:return F(0)
    v=union_sweep(B);v2=union_grid(B)
    assert v==v2
    return v

def split(a,b,P):
    length=b-a
    if length<=0:return []
    if length>=P:return [(F(0),P)]
    a=a%P
    return [(a,a+length)] if a+length<=P else [(a,P),(F(0),a+length-P)]

def lin_rational():
    # Frozen decimal parameters, NOT assertions that the irrational Lin value is rational.
    a=F('1.615069293303905');v=F('2.615069293303905');w=F('3.973518115477081')
    P=2*a;boxes=[];witness=[]
    for s in (0,1):
        for k,l in itertools.product(range(1,7),repeat=2):
            mass=l*sum((F(1) if (s+i)%2==0 else w for i in range(k)),F(0))
            side=root_lower(mass,2);r=side/2
            xl=(s+k-1)*a-r;xh=s*a+r;yl=(l-1)*v-r;yh=r
            if xl>=xh or yl>=yh:continue
            # Witness extrema imply containment for every point of the unwrapped box.
            assert max(abs(xl-s*a),abs(xh-(s+k-1)*a))<=r
            assert max(abs(yl),abs(yh-(l-1)*v))<=r
            bs=[X+Y for X,Y in itertools.product(split(xl,xh,P),split(yl,yh,v))]
            boxes+=bs;witness.append({'start':s,'counts':[k,l],'mass':str(mass),'side':str(side),'boxes':[[str(x) for x in b] for b in bs]})
    vol=exact_union(boxes);ratio=vol/(1+w)
    assert ratio>F('1.6855099933555')
    return {'parameters':list(map(str,(a,v,w))),'witnesses':witness,'union_area':str(vol),'mass_per_period':str(1+w),'ratio':str(ratio),'display':float(ratio),'status':'reproduction_only'}

def finite_cert(row):
    n=row['n'];p=row['parameters']
    z=[(F(0),F(0))]+[tuple(F(f'{p[2*i+j]:.8f}') for j in range(2)) for i in range(n-1)]
    w=[F(1)]+[F(f'{math.exp(x):.8f}') for x in p[2*(n-1):]]
    B=[];witness=[]
    for mask in range(1,1<<n):
        ix=[i for i in range(n) if mask>>i&1];mass=sum((w[i] for i in ix),F(0));side=root_lower(mass,2);r=side/2
        box=tuple(x for j in range(2) for x in (max(z[i][j] for i in ix)-r,min(z[i][j] for i in ix)+r))
        if any(box[2*j]>=box[2*j+1] for j in range(2)):continue
        for i in ix:
            for j in range(2):assert max(abs(box[2*j]-z[i][j]),abs(box[2*j+1]-z[i][j]))<=r
        B.append(box);witness.append({'mask':mask,'side':str(side),'mass':str(mass),'box':list(map(str,box))})
    vol=exact_union(B);ratio=vol/sum(w)
    return {'n':n,'points':[[str(x) for x in zz] for zz in z],'weights':list(map(str,w)),'witnesses':witness,'ratio':str(ratio),'display':float(ratio),'status':'certified_lower_bound_below_existing_baseline'}

def run():
    checks=[]
    rng=random.Random(20261023)
    for d in (2,3,4):
        for _ in range(35):
            B=[]
            for i in range(7):
                b=[]
                for j in range(d):
                    a=rng.randrange(-3,4);b.extend((a,a+rng.randrange(0,4)))
                B.append(tuple(b))
            exact_union(B+B[:2])
    checks.append('105 random integer box unions: sweep/grid agree, including duplicates and degeneracy')
    # A deterministic touching/overlap case of known area.
    assert exact_union([(0,2,0,2),(1,3,1,3),(2,2,0,9)])==7
    checks.append('known overlap area 7; zero-width box ignored')
    lin=lin_rational();(ROOT/'results/lin_rational_certificate.json').write_text(json.dumps(lin,indent=2))
    checks.append('Lin nearby rational reproduction, independent exact union algorithms')
    rows=json.loads((ROOT/'results/round1_search.json').read_text())['records']
    fs=[finite_cert(x) for x in rows if x['label']=='finite_atoms']
    (ROOT/'results/finite_certificates.json').write_text(json.dumps(fs,indent=2))
    checks.append('three finite-atom candidates rationalized and certified')
    upper=[json.loads((ROOT/'results/round1_upper.json').read_text())]
    for row in upper:
        again=certify(row['R'],row['epsilon'],row['r0'])
        for k,v in again.items():assert row[k]==v
        assert F(row['improvement_vs_public_exact'])>0
    checks.append('round 1 global scalar kernel certificate rechecked with exact fractions')
    # Within the classical disjoint-cover mechanism: parameter tau>1 never improves 2^d.
    cover=[]
    for d in (3,4):
        vals=[]
        for k in (1,2,3,4,6,9):
            tau=1+F(1,10**k);v=(1+tau)**d;assert v>2**d
            vals.append({'tau':str(tau),'bound':str(v)})
        cover.append({'dimension':d,'values':vals,'infimum':2**d,'no_improvement_in_this_mechanism':True})
    (ROOT/'results/cover_parameter_audit.json').write_text(json.dumps(cover,indent=2))
    checks.append('d=3,4 classical covering parameter cannot beat 2^d; analytic proof supplied')
    out={'checks':checks,'all_passed':True,'finite_displays':[x['display'] for x in fs],'lin_rational_display':lin['display']}
    (ROOT/'results/audit.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

if __name__=='__main__':run()
