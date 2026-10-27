#!/usr/bin/env python3
"""Independent real-variable ALL-CONE density and local-series integral audit.

Unlike the author producer, no projected-node groups or z=-it accumulation
are used. All 1152 real and imaginary numerators are added over one fixed
common denominator in u=t^2 before any cancellation. Partial fractions use
local Taylor quotient recurrences rather than SymPy apart/residue.
"""
from fractions import Fraction as Q
from collections import Counter
from itertools import permutations,product,combinations
from math import comb,isqrt,prod
from pathlib import Path
import datetime,hashlib,json,resource,sys,time
from flint import fmpq_poly as Poly,fmpq
import flint
import sympy as sp

D=Path(__file__).resolve().parent
A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube8-terra-scalar/cube-dim-8-bounds/upper/round004/block_green_anchors')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def qp(p):return [Q(str(c)) for c in p]
def lin(b):return Poly([fmpq(str(b)),1])
def poly_product(polys):
    out=Poly([1])
    for p in polys:out*=p
    return out
def determinant(rows):
    a=[[Q(c) for c in row] for row in rows];out=Q(1)
    for k in range(len(a)):
        p=next(j for j in range(k,len(a)) if a[j][k])
        if p!=k:a[k],a[p]=a[p],a[k];out=-out
        out*=a[k][k]
        for j in range(k+1,len(a)):
            ratio=a[j][k]/a[k][k]
            for l in range(k,len(a)):a[j][l]-=ratio*a[k][l]
    return out
def dot(v,w):return sum(a*b for a,b in zip(v,w))

def build_cones():
    body=list(product((1,),(-1,1),(-1,1),(-1,1)))
    rays=[(1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1),(2,1,1,0),(3,1,1,1)]
    ids=[(0,1,4,5),(4,1,2,3),(4,1,3,5)]
    # Explicit inverse rows, checked against each basis column matrix.
    inverses=[[(1,-1,-1,-1),(0,1,-1,0),(0,0,1,-1),(0,0,0,1)],
              [(1,-1,0,0),(0,1,-1,0),(-1,1,1,-1),(0,0,0,1)],
              [(0,0,1,-1),(0,1,-1,0),(Q(-1,2),Q(1,2),Q(1,2),Q(1,2)),(Q(1,2),Q(-1,2),Q(-1,2),Q(1,2))]]
    base=[]
    for inds,inverse,expected_det in zip(ids,inverses,[1,1,2]):
        cone=tuple(rays[j] for j in inds)
        assert abs(determinant(cone))==expected_det
        assert [[dot(row,col) for col in cone] for row in inverse]==[[int(i==j) for j in range(4)] for i in range(4)]
        base.append((cone,expected_det))
    cones=[]
    # Act by sending old coordinate j to position permutation[j]. This is
    # an independently chosen action convention; the entire group is used.
    for perm in permutations(range(4)):
        for signs in product((-1,1),repeat=4):
            for cone,det in base:
                transformed=[]
                for ray in cone:
                    v=[0]*4
                    for j in range(4):v[perm[j]]=signs[j]*ray[j]
                    transformed.append(tuple(v))
                cones.append((tuple(transformed),det))
    assert len(cones)==1152 and len({tuple(sorted(c)) for c,d in cones})==1152
    allrays={r for c,d in cones for r in c};assert len(allrays)==240
    rate={r:sum(abs(dot(v,r)) for v in body) for r in allrays}
    assert [rate[r] for r in rays]==[8,8,12,12,16,24]
    for cone,det in cones:
        for v in body:
            signs=[dot(v,r) for r in cone]
            assert min(signs)>=0 or max(signs)<=0
    raw_weight=sum((Q(det,16*prod(rate[r] for r in cone)) for cone,det in cones),Q(0))
    assert raw_weight==Q(11,3072)
    return cones,rate,raw_weight

def density(point,cones,rate):
    raw=[];common_powers=Counter();zero_cones=0;repeated_cones=0
    for cone,det in cones:
        bs=[Q(dot(point,r),rate[r]) for r in cone]
        poles=Counter(b*b for b in bs)
        for pole,m in poles.items():common_powers[pole]=max(common_powers[pole],m)
        zero_cones+=int(Q(0) in poles)
        repeated_cones+=int(len(set(bs))<4)
        e1=sum(bs,Q(0));e2=sum((x*y for x,y in combinations(bs,2)),Q(0))
        e3=sum((prod(triple) for triple in combinations(bs,3)),Q(0));e4=prod(bs)
        real=Poly([fmpq(str(e4)),fmpq(str(-e2)),1])
        imaginary_over_t=Poly([fmpq(str(-e3)),fmpq(str(e1))])
        raw.append((poles,Q(det,16*prod(rate[r] for r in cone)),real,imaginary_over_t))
    H=poly_product([lin(pole)**m for pole,m in common_powers.items()])
    numerator=Poly([0]);imaginary=Poly([0])
    for poles,weight,real,im in raw:
        rawden=poly_product([lin(pole)**m for pole,m in poles.items()])
        multiplier,remainder=divmod(H,rawden);assert not remainder
        numerator+=fmpq(str(weight))*real*multiplier
        imaginary+=fmpq(str(weight))*im*multiplier
    assert not imaginary
    # First cancellation occurs here, after all 1152 terms were summed.
    common=numerator.gcd(H);N=numerator//common;H=H//common
    assert H[0]!=0 and N[0]!=0
    # Normalize denominator leading coefficient to 1 for stable recording.
    lead=H.leading_coefficient();N=N/lead;H=H/lead
    multiplicities={};remaining=H
    for pole in sorted(common_powers):
        m=0
        while remaining.degree()>0:
            quotient,remainder=divmod(remaining,lin(pole))
            if remainder:break
            m+=1;remaining=quotient
        if m:multiplicities[pole]=m
    assert remaining.degree()==0 and Q(0) not in multiplicities
    assert H.degree()-N.degree()==2
    leading=Q(str(N.leading_coefficient()/H.leading_coefficient()))
    assert leading==Q(11,3072)
    return N,H,multiplicities,{'raw_common_denominator_degree':sum(common_powers.values()),
      'raw_zero_pole_multiplicity':common_powers.get(Q(0),0),
      'cones_with_zero_projection':zero_cones,'cones_with_repeated_signed_projection':repeated_cones,
      'complete_imaginary_numerator_zero':True,'reduced_numerator_degree_u':N.degree(),
      'reduced_denominator_degree_u':H.degree(),'t0_pi_minus4':str(Q(str(N[0]/H[0]))),
      'leading_t_minus4_coefficient_pi_minus4':str(leading),
      'N_coefficients_ascending_u':list(map(str,qp(N))),'D_coefficients_ascending_u':list(map(str,qp(H))),
      'pole_squared_multiplicities':{str(p):m for p,m in multiplicities.items()}}

def shift_coefficients(poly,root,count):
    # Taylor coefficients of poly(root+v), via the binomial theorem.
    coeffs=qp(poly)
    return [sum((coeffs[j]*comb(j,k)*root**(j-k) for j in range(k,len(coeffs))),Q(0)) for k in range(count)]

def integrate_rational(N,H,poles):
    assert N.degree()<H.degree()
    reconstructed=Poly([0]);coefficient=Q(0);terms=[]
    for squared,m in sorted(poles.items()):
        assert squared>0
        sn,sd=isqrt(squared.numerator),isqrt(squared.denominator)
        assert sn*sn==squared.numerator and sd*sd==squared.denominator
        a=Q(sn,sd)
        other,remainder=divmod(H,lin(squared)**m);assert not remainder
        numerator=shift_coefficients(N,-squared,m)
        denominator=shift_coefficients(other,-squared,m);assert denominator[0]
        expansion=[]
        for j in range(m):
            expansion.append((numerator[j]-sum((denominator[k]*expansion[j-k] for k in range(1,j+1)),Q(0)))/denominator[0])
        for j,c in enumerate(expansion):
            k=m-j
            if not c:continue
            termden=lin(squared)**k
            multiplier,remainder=divmod(H,termden);assert not remainder
            reconstructed+=fmpq(str(c))*multiplier
            coefficient+=c*Q(comb(2*k-2,k-1),2*4**(k-1))/a**(2*k-1)
            terms.append({'a':str(a),'k':k,'coefficient':str(c)})
    assert reconstructed==N
    return coefficient,terms

def parse_author_density(text):
    t=sp.Symbol('t');expr=sp.sympify(text,locals={'t':t})
    N,H=map(lambda x:sp.Poly(x,t),sp.fraction(expr))
    assert all(e[0]%2==0 for e,c in N.terms()+H.terms())
    def convert(p):
        out=[Q(0)]*(p.degree()//2+1)
        for (power,),c in p.terms():out[power//2]=Q(c)
        return Poly([fmpq(str(c)) for c in out])
    return convert(N),convert(H)

def main():
    author=json.loads((A/'result.json').read_text());manifest=json.loads((A/'manifest.json').read_text())
    for name,digest in manifest['files'].items():assert sha(A/name)==digest,name
    for receipt_name in ['receipt.json','receipt_interrupted.json']:
        receipt=json.loads((A/receipt_name).read_text())
        for path,digest in receipt['source_sha256'].items():
            p=Path(path)
            # The interrupted receipt records its then-current run.py path;
            # the exact old bytes are now intentionally run_interrupted.py.
            if receipt_name=='receipt_interrupted.json' and p==A/'run.py':p=A/'run_interrupted.py'
            assert sha(p)==digest,str(p)
    cones,rate,weight=build_cones()
    forms={};rows={};four_dimensional={}
    for name,point in [('A',(1,1,1,0)),('B',(1,1,0,0))]:
        N,H,poles,row=density(point,cones,rate)
        savedN,savedH=parse_author_density(author['p4_pi_minus4'][name])
        assert N*savedH==H*savedN
        coef,pf=integrate_rational(N,H,poles)
        assert coef=={'A':Q(17,640),'B':Q(1,24)}[name]
        forms[name]=(N,H,poles);rows[name]=row
        four_dimensional[name]={'G4_pi_minus3':str(coef),'density_full_function_identity_verified':True}
    targets={'AA':Q(10495703317,2149908480000),'AB':Q(15547486331,1920360960000),'BB':Q(143,10584)}
    integrals={}
    for label in ['AA','AB','BB']:
        (N,H,poles),(n,h,otherpoles)=[forms[c] for c in label]
        NN,HH=N*n,H*h
        merged=Counter(poles);merged.update(otherpoles)
        coef,pf=integrate_rational(NN,HH,merged)
        assert coef==targets[label]
        source=next(row for row in author['anchors'] if row['point']==label)
        assert Q(source['G8_pi_minus7_coefficient'])==coef
        # Verify the author's PF identity separately by cross multiplication.
        sourceN=Poly([0])
        for item in source['partial_fractions']:
            pole=Q(item['a'])**2;k=item['k'];c=Q(item['coefficient'])
            multiplier,remainder=divmod(HH,lin(pole)**k);assert not remainder
            sourceN+=fmpq(str(c))*multiplier
        assert sourceN==NN
        integrals[label]={'G8_pi_minus7':str(coef),'partial_fraction_terms_independent_local_Taylor':pf,
          'full_partial_fraction_identity':True,'source_partial_fraction_identity':True,
          't0_integrand_pi_minus8':str(Q(str(NN[0]/HH[0]))),
          'tail_bound_pi_minus8_coefficient':str(weight*weight),
          'tail_integral_T_power':-7,'tail_integral_coefficient_pi_minus8':str(weight*weight/7)}
    result={'status':'PASS','scope':'three exact point anchors only; no g_min, continuous I certificate or new C8 bound',
      'independent_density_method':'all 1152 real/imaginary numerators over one common u=t^2 denominator, first gcd only after complete summation',
      'cone_count':len(cones),'ray_count':len(rate),'density_prefactor':'1/(16*pi^4)',
      'shared_time_model':'p8_t(x,y)=p4_t(x)*p4_t(y), integrated over one t',
      'raw_absolute_density_tail_coefficient_pi_minus4':str(weight),
      'exact_global_peak':{'p4_at_origin_coefficient_pi_minus4_t_minus4':str(weight),
       'p8_at_origin_coefficient_pi_minus8_t_minus8':str(weight*weight),
       'global_uniform_long_tail_coefficient_pi_minus8_T_minus7':str(weight*weight/7),
       'improvement_factor_over_basis_C':str(Q(1,256)/(weight*weight)),
       'proof':'At x=0 every raw cone term is positive weight*t^-4; positive Fourier inversion bounds abs(p4_t(x)) by p4_t(0). Shared-time product gives the p8 peak.'},
      'densities':rows,'known_four_dimensional_integral_crosschecks':four_dimensional,
      'eight_dimensional_integrals':integrals,
      'author_quadrature_repeated':False,
      'author_failed_run_preserved':{'status':'KeyboardInterrupt; exit 1','wall_seconds':json.loads((A/'receipt_interrupted.json').read_text())['wall_seconds'],
       'script':'run_interrupted.py','receipt':'receipt_interrupted.json'},
      'author_files':{p.name:sha(p) for p in sorted(A.iterdir()) if p.is_file()}}
    (D/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['densities','eight_dimensional_integrals','author_files']},indent=2))
    print(json.dumps({'point_coefficients':{k:v['G8_pi_minus7'] for k,v in integrals.items()},'density_metadata':{k:{j:v for j,v in row.items() if not j.endswith('ascending_u')} for k,row in rows.items()}},indent=2))

if __name__=='__main__':
    start=datetime.datetime.now(datetime.timezone.utc).isoformat();wall=time.perf_counter();cpu=time.process_time();code=0;error=None
    try:main()
    except BaseException as exc:code=1;error=repr(exc);raise
    finally:
        receipt={'status':'PASS' if code==0 else 'FAIL','exit_code':code,'error':error,
         'start_utc':start,'finish_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'wall_seconds':time.perf_counter()-wall,'cpu_seconds':time.process_time()-cpu,
         'peak_rss_native_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
         'python':sys.version,'executable':sys.executable,'flint':flint.__version__,'sympy':sp.__version__,
         'checker_sha256':sha(Path(__file__)),
         'dependencies':'Fraction; flint exact univariate rational polynomials; SymPy only parses author polynomial strings'}
        (D/'actualreceipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
