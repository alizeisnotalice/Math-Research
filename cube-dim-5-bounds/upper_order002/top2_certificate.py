#!/usr/bin/env python3
"""Single point, one registered 256-bin top-two order certificate."""
from pathlib import Path
from itertools import product,combinations
from collections import Counter
from flint import arb,ctx
import sympy as sp
import hashlib,json,sys,time

HERE=Path(__file__).resolve().parent

def stats(z,terms=64):
    beta=arb(3)/4
    assert z>=arb(1)/4
    h=(arb(1)/2).sqrt()
    starts=[beta.gamma()*h/arb.pi()*z**(-beta),
            (2*beta).gamma()/(2*arb.pi())*z**(-2*beta),
            (3*beta).gamma()*h/(6*arb.pi())*z**(-3*beta)]
    tail=arb(0);M=arb(0);Itail=arb(0)
    for residue,term in enumerate(starts,1):
        st=arb(0);sm=arb(0);sn=arb(0)
        for j in range(terms):
            k=residue+4*j;st+=term
            sm+=term*(beta*k)/(1-beta*k)*z
            sn+=term*(beta*k)/(1+beta*k)/z
            term=-term*(beta*k)*(beta*k+1)*(beta*k+2)*z**(-3)/((k+1)*(k+2)*(k+3)*(k+4))
        k=residue+4*terms
        assert z**(-3)/k<1
        tail+=st+arb(0).union(term)
        M+=sm+arb(0).union(term*(beta*k)/(1-beta*k)*z)
        Itail+=sn+arb(0).union(term*(beta*k)/(1+beta*k)/z)
    F=1-tail;n0=(1/beta).gamma()/beta;m=M/F;a=F/(n0-Itail)
    assert F>0 and F<1 and m>a and m<z
    return F,m,a

def det(c,z,w,r2):
    return (121-r2)*c*c+(55+r2)*c*(z+w)+(25-r2)*z*w

def coeffs(key):
    r2,B,D,C=key
    return (484-4*r2+11*(B+D)+2*C,
            220+4*r2-11*B+5*D-2*C,
            220+4*r2+5*B-11*D-2*C,
            100-4*r2-5*(B+D)+2*C)

def quadratic(a,z,w,r2,co):
    return (co[0]*a*a+co[1]*a*z+co[2]*a*w+co[3]*z*w)/(16*a*det(a,z,w,r2))

def main():
    start=time.monotonic();ctx.prec=256
    directions=[(1,)+t for t in product((-1,1),repeat=4)]
    classes=Counter();representatives={}
    for v,w in combinations(directions,2):
        r=sum(x*y for x,y in zip(v,w));b=sum(v[:4]);d=sum(w[:4])
        if b*b<d*d:v,w=w,v;b,d=d,b
        key=(r*r,b*b,d*d,r*b*d);classes[key]+=1;representatives[key]=(v,w)
    assert len(classes)==11 and sum(classes.values())==120
    c,z,w=sp.symbols('c z w',positive=True);x=sp.Matrix([1,1,1,1,0])
    # Exact identity checks for one integer pair in every class.
    for key,(vv,ww) in representatives.items():
        v=sp.Matrix(vv);v2=sp.Matrix(ww);co=coeffs(key)
        assert all(q>=0 for q in co)
        B=c*(16*sp.eye(5)-v*v.T-v2*v2.T)+z*v*v.T+w*v2*v2.T
        assert sp.simplify(B.det()-(16*c)**3*det(c,z,w,key[0]))==0
        assert sp.simplify((x.T*B.inv()*x)[0]-quadratic(c,z,w,key[0],co))==0
    L,H,t=sp.symbols('L H t',positive=True)
    assert sp.simplify(2*sp.integrate((H-t)*t**14,(t,L,H))-2*(H*(H**15-L**15)/15-(H**16-L**16)/16))==0
    beta=arb(3)/4;n0=(1/beta).gamma()/beta
    kappa=3*(arb(7)/4).gamma()/(arb(2)**(arb(7)/2)*arb.pi()**(arb(5)/2))
    endpoints=[arb(1)/4]
    for _ in range(256):endpoints.append(endpoints[-1]*17/16)
    moments=[stats(t) for t in endpoints]
    deltaF=[b[0]-a[0] for a,b in zip(moments,moments[1:])]
    weights_min=[(b[0]**15-a[0]**15)/15 for a,b in zip(moments,moments[1:])]
    diagonal=[2*(b[0]*(b[0]**15-a[0]**15)/15-(b[0]**16-a[0]**16)/16) for a,b in zip(moments,moments[1:])]
    assert all(w>0 for w in deltaF+weights_min+diagonal)
    keys=sorted(classes);coefficient=[coeffs(k) for k in keys]
    lows=[arb(0) for k in keys];highs=[arb(0) for k in keys]
    weight_sum=arb(0);completed_rows=0
    for i in range(256):
        if time.monotonic()-start>120:
            checkpoint={'status':'INCOMPLETE_CHECKPOINT','completed_rows':completed_rows,'registered_rows':256,'seconds':time.monotonic()-start,'partial_lower_groups':[str(t) for t in lows],'partial_upper_groups':[str(t) for t in highs],'scope':'No complete formula enclosure; source can replay registered batch.'}
            (HERE/'checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
            raise TimeoutError('Registered120s checkpoint reached')
        zlo,zhi=endpoints[i:i+2]
        for j in range(i+1):
            wlo,whi=endpoints[j:j+2]
            F_lo,mlo,alo=moments[j];F_hi,mhi,ahi=moments[j+1]
            probability=diagonal[j] if i==j else deltaF[i]*weights_min[j]
            weight_sum+=probability if i==j else 2*probability
            common={}
            for r2 in (1,9):
                detlow=(16*mhi)**(-arb(3)/2)*det(mhi,zhi,whi,r2)**(-arb(1)/2)
                dethigh=(16*mlo)**(-arb(3)/2)*det(mlo,zlo,wlo,r2)**(-arb(1)/2)
                common[r2]=(detlow,dethigh)
            for k,key in enumerate(keys):
                r2=key[0];co=coefficient[k];dl,dh=common[r2]
                qlo=quadratic(alo,zlo,wlo,r2,co);qhi=quadratic(ahi,zhi,whi,r2,co)
                low=dl*qlo**(-arb(7)/4);high=dh*qhi**(-arb(7)/4)
                if i!=j:
                    qlo2=quadratic(alo,wlo,zlo,r2,co);qhi2=quadratic(ahi,whi,zhi,r2,co)
                    low+=dl*qlo2**(-arb(7)/4);high+=dh*qhi2**(-arb(7)/4)
                assert high>low and low>0
                lows[k]+=probability*low;highs[k]+=probability*high
        completed_rows=i+1
    FL=moments[0][0];FU,mU,aU=moments[-1];ell=endpoints[0];U=endpoints[-1]
    exact_mass=FU**16-16*FU*FL**15+15*FL**16
    assert (120*weight_sum-exact_mass).contains(0)
    # Uniform kernel bound: C0*m(t)^(-3/2)*t^(5/4)*max^(-1/2).
    C0=arb(20)**(arb(7)/4)/(4*arb(16)**(arb(3)/2))
    low_omission_per_pair=2*C0*n0**(arb(3)/2)*ell**beta*FL**(arb(27)/2)/(arb(27)/2)
    A=arb(0)
    for j,probability in enumerate(weights_min):
        A+=probability*moments[j][1]**(-arb(3)/2)*endpoints[j+1]**(arb(5)/4)
    one_above_U=2*C0*U**(-arb(1)/2)*(1-FU)*A
    Ctail=1/(1-(-arb(1)).exp())
    both_above_U=C0*mU**(-arb(3)/2)*(U**beta*(1-FU)**2+Ctail**2*U**(-beta))
    low_omission=120*kappa*low_omission_per_pair
    high_omission=120*kappa*(one_above_U+both_above_U)
    output=[];lower=arb(0);upper_finite=arb(0)
    for key,lo,hi in zip(keys,lows,highs):
        count=classes[key];a=kappa*count*lo;b=kappa*count*hi
        lower+=a;upper_finite+=b
        output.append({'invariants':list(key),'count':count,'lower':str(a),'finite_upper':str(b)})
    upper=upper_finite+low_omission+high_omission
    result={'status':'TOP2_ORDER_SINGLE_POINT_CERTIFIED_INTERVAL','point':[1,1,1,1,0],'entropy_penalty':0,'unordered_pair_count':120,'class_count':11,
      'symbolic_checks':{'all11_determinants_and_quadratics':True,'quadratic_positive_polynomial_coefficients':True,'diagonal_weight':True,'partition_mass_closure':True},
      'partition':{'left':'1/4','ratio':'17/16','bins':256,'right':str(U),'included_top2_probability':str(exact_mass),'all_weights_positive':True},
      'class_contributions':output,'Green_lower_Arb':str(lower),'Green_formula_upper_Arb':str(upper),'Green_low_omission_upper':str(low_omission),'Green_high_omission_upper':str(high_omission),
      'strictly_improves_order001_upper':lower>arb('0.0000431012868806624'),'certifies_1_over10000':lower>arb(1)/10000,'formula_strictly_below_1_over10000':upper<arb(1)/10000,
      'scope':'Full top2 conditional Jensen formula at11110; upper endpoint does not bound actualG. No top3 or spatial grid.',
      'precision_bits':ctx.prec,'seconds':time.monotonic()-start,'python':sys.executable,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
