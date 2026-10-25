"""Independent exact g audit: Leibniz, direct pairs, closed coefficient moments."""
from fractions import Fraction as F
from itertools import combinations,permutations,product
from math import comb,factorial,prod
from pathlib import Path
import hashlib,json,sys,time
import mpmath as mp

ROOT=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds')
BETA=F(3,4)
V=[(1,)+x for x in product((-1,1),repeat=4)]
PERMS=[(p,(-1)**sum(p[i]>p[j] for i in range(5) for j in range(i+1,5)))
       for p in permutations(range(5))]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def determinant(rows):
    return sum(sign*prod(rows[i][p[i]] for i in range(5)) for p,sign in PERMS)

def rising(beta,n):return prod((beta+j for j in range(n)),start=F(1))

def weight_moment(k,ell):
    # Coefficient reciprocal expansion, with k variables squared and ell linear.
    n=2*k+ell
    return sum((F(comb(k,j)*factorial(ell+k+j),16**(ell+k+j))
                *BETA**(ell+k+j)*(1-BETA)**(k-j) for j in range(k+1)),F(0))/rising(BETA,n)

def geometry(size,deadline):
    items=[];total=0;enumerated=0
    for subset in combinations(range(16),size):
        rows=[V[j] for j in subset]
        if size==4:rows.append((1,1,1,1,1))
        w=determinant(rows)**2
        enumerated+=1;total+=w
        if w:items.append((sum(1<<j for j in subset),w))
    overlaps=[0]*(size+1)
    for i,(mask,w) in enumerate(items):
        overlaps[size]+=w*w
        for j in range(i):
            mask2,w2=items[j]
            overlaps[(mask&mask2).bit_count()]+=2*w*w2
        if time.monotonic()>deadline:raise TimeoutError('120 second finite batch limit')
    assert sum(overlaps)==total**2
    moments=[weight_moment(k,2*size-2*k) for k in range(size+1)]
    mean=total*weight_moment(0,size)
    second=sum((a*b for a,b in zip(overlaps,moments)),F(0))
    return {'subsets':enumerated,'nonzero':len(items),'coefficient_sum':total,
            'ordered_pair_overlap_weights':overlaps,
            'overlap_moments':[str(x) for x in moments],
            'mean':str(mean),'second_moment':str(second)}

def atan_bounds(q,n=20):
    lower=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(n)),F(0))
    return lower,lower+F(1,(2*n+1)*q**(2*n+1))

def main():
    began=time.monotonic();deadline=began+120
    sources=[AUTHOR/d/name for d in ('astra_upper001','astra_upper002')
             for name in ('registration.md','report.md','verify.py','verification.json')]
    hashes={str(p):sha(p) for p in sources}
    old=json.loads((AUTHOR/'astra_upper001/verification.json').read_text())
    target=json.loads((AUTHOR/'astra_upper002/verification.json').read_text())
    # Author code is data only: never imported or executed.
    for path,expected in target['source_sha256'].items():
        assert sha(Path(path))==expected
    frame=[[sum(v[i]*v[j] for v in V) for j in range(5)] for i in range(5)]
    assert frame==[[16*int(i==j) for j in range(5)] for i in range(5)]
    d=geometry(5,deadline);b=geometry(4,deadline)
    assert d['subsets']==4368 and d['nonzero']==3008 and d['coefficient_sum']==16**5
    assert b['subsets']==1820 and b['nonzero']==940 and b['coefficient_sum']==5*16**4
    for receipt,key in ((d,'finite_determinant_minorant'),):
        assert receipt['ordered_pair_overlap_weights']==old[key]['ordered_pair_overlap_weight_sums']
        assert receipt['overlap_moments']==old[key]['overlap_moments']
    assert b['ordered_pair_overlap_weights']==target['cofactor_receipt']['ordered_pair_overlap_weights']
    assert b['overlap_moments']==target['cofactor_receipt']['overlap_weight_moments']
    m1,m2,B1,B2=map(F,(d['mean'],d['second_moment'],b['mean'],b['second_moment']))
    assert (m1,m2,B1,B2)==(F(648,1463),F(126404864,474619145),F(216,77),F(1517312,149017))
    assert m2/m1==F(15800608,26277615)
    a5=atan_bounds(5);a239=atan_bounds(239)
    pil,pih=16*a5[0]-4*a239[1],16*a5[1]-4*a239[0]
    assert F(314159,100000)<pil<pih<F(355,113)
    # g3^8 = 3^8 m1^20/(2^60 pi^20 m2^5 B2^7).
    ratio_to_rational_target_eighth=F(3)**8*m1**20*27000**8/(2**60*F(355,113)**20*m2**5*B2**7)
    assert ratio_to_rational_target_eighth>1
    # Ceiling^4 = 3^4 m1^5/(2^30 pi^10 B1^7).
    ceiling_ratio_fourth=F(3)**4*m1**5*10000**4/(2**30*F(314159,100000)**10*B1**7)
    assert ceiling_ratio_fourth<1
    mp.mp.dps=60
    def dec(x):return mp.mpf(x.numerator)/x.denominator
    c=3/(mp.power(2,mp.mpf(15)/2)*mp.pi**(mp.mpf(5)/2))
    g=c*dec(m1)**(mp.mpf(5)/2)/dec(m2)**(mp.mpf(5)/8)/dec(B2)**(mp.mpf(7)/8)
    kappa=(4*mp.pi)**(-mp.mpf(5)/2)*(mp.mpf(3)/4)*mp.power(4,mp.mpf(7)/4)*mp.gamma(mp.mpf(7)/4)
    c_from_gamma=kappa/(16*mp.gamma(mp.mpf(7)/4))
    assert mp.almosteq(c,c_from_gamma)
    out={'status':'independent_exact_global_g_1_over_27000_audit_pass',
         'geometry_determinant':d,'geometry_adjugate_corner':b,'frame':frame,
         'exact_inputs':{'Ed':str(m1),'Ed2':str(m2),'Eb_corner':str(B1),'Eb2_corner':str(B2)},
         'strict_rational_comparison':{
             'g3_target':'1/27000','pi_upper':'355/113','pi_lower':'314159/100000',
             'target_ratio_eighth_exact':str(ratio_to_rational_target_eighth),
             'target_ratio_eighth_gt_one':True,
             'holder_ceiling_ratio_fourth_exact':str(ceiling_ratio_fourth),
             'holder_ceiling_lt_1_over_10000':True,
             'pi_bounds_verified_by':'20-term alternating Machin arctan rational intervals'},
         'mpmath_display_only':{'g3':str(g),'c':str(c),'c_from_time_Gamma_cancellation':str(c_from_gamma)},
         'limitations':['Not a new C5 bound or an I certificate','No author arithmetic engine imported',
                         'Fraction rational certificate independent of shared Arb; mpmath decimals not certified'],
         'author_input_sha256':hashes,'source_sha256':sha(Path(__file__)),
         'registration_sha256':sha(ROOT/'registration.md'),
         'runtime':{'python':sys.executable,'seconds':time.monotonic()-began}}
    (ROOT/'independent_check.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','exact_inputs','mpmath_display_only','runtime','source_sha256')},indent=2))

if __name__=='__main__':main()
