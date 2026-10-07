"""Fixed symmetric jump absorption identities, exact rational + numerical energy.

No actual cube/FIRST input is constructed. Numeric quadrature is not interval
certified; omitted tails have separate analytic finite-matrix bounds.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np

BASE = Path(__file__).parent
REG = BASE/'nonlocal_abel_guard_registration_20261007.json'
COUNT = 0

def check(v, name):
    global COUNT
    COUNT += 1
    if not v:
        raise AssertionError(name)

def mv(a,v):
    return [sum((x*y for x,y in zip(row,v)),F()) for row in a]

def add(v,w):
    return [x+y for x,y in zip(v,w)]

def sub(v,w):
    return [x-y for x,y in zip(v,w)]

def scale(v,c):
    return [c*x for x in v]

def resolvent(a,s):
    m = len(a)
    rows = [[a[i][j]+(s if i==j else 0) for j in range(m)]
            +[F(int(i==j)) for j in range(m)] for i in range(m)]
    for j in range(m):
        pivot = next(i for i in range(j,m) if rows[i][j])
        rows[j],rows[pivot] = rows[pivot],rows[j]
        p = rows[j][j]
        rows[j] = [v/p for v in rows[j]]
        for i in range(m):
            if i != j:
                c = rows[i][j]
                rows[i] = [x-c*y for x,y in zip(rows[i],rows[j])]
    return [[s*v for v in row[m:]] for row in rows]

def powers(a,v,N):
    out = [v]
    for _ in range(N):
        out.append(mv(a,out[-1]))
    return out

def receipt(v):
    # Length-prefixed binary representation avoids decimal serialization limits.
    h = hashlib.sha256()
    for x in v:
        for k in (x.numerator,x.denominator):
            b = abs(k).to_bytes(max(1,(abs(k).bit_length()+7)//8),'big')
            h.update(bytes([k<0])); h.update(len(b).to_bytes(8,'big')); h.update(b)
    return {'length':len(v),'sha256':h.hexdigest()}

def run(m,d,qorder):
    before = COUNT
    J = [[F(1+(i+j)%3,8*m) if i!=j else F() for j in range(m)] for i in range(m)]
    degree = list(map(sum,J))
    L = [[degree[i] if i==j else -J[i][j] for j in range(m)] for i in range(m)]
    u = [F(1+i%3,64*m) if i<d else F() for i in range(m)]
    Lu = mv(L,u)
    nu = [1+Lu[i] if i<d else F() for i in range(m)]
    mu = sub(nu,Lu)
    V = [1/u[i] for i in range(d)]
    K = [[L[i][j]+(V[i] if i==j else 0) for j in range(d)] for i in range(d)]
    ext = lambda z: z+[F()]*(m-d)
    def T(z):
        return [V[i]*z[i] if i<d else sum((J[i][j]*z[j] for j in range(d)),F()) for i in range(m)]
    check(all(x>=0 for x in nu),'positive original source')
    check(all(0<=x<=1 for x in mu),'whole-space cap')
    check(sum(nu)==sum(mu),'full mass conservation')
    check(T(u[:d])==mu,'inside plus outside absorption of odometer')
    check(mv(K,u[:d])==nu[:d],'killed inverse source')
    check(sum(mu[d:])>0,'genuine exterior jump flux')
    exact=[]
    for s in (F(1,4),F(1),F(4)):
        A,B = resolvent(L,s),resolvent(K,s)
        Anu,Bnu = powers(A,nu,9),powers(B,nu[:d],9)
        Bu = powers(B,u[:d],10)
        z0 = sub(u[:d],Bu[1])
        zj = powers(B,z0,9)
        rho = [T(z) for z in zj]
        for j,z in enumerate(zj):
            check(mv(L,ext(z))==sub(ext(mv(K,z)),T(z)),'zero-extension leakage identity')
            check(all(0<=x<=y for x,y in zip(rho[j],mu)),'positive layer cap')
            check(all((j+1)*x*x<=y*y for x,y in zip(rho[j],mu)),'Poisson peak layer envelope')
        check(add([sum(col) for col in zip(*rho[:9])],T(Bu[9]))==mu,'exact infinite-sum remainder after nine layers')
        for N in (1,2,4,8):
            allpowers = [powers(A,rho[j],N-j+1) for j in range(N)]
            rhs = [sum(items) for items in zip(*(row[N-j] for j,row in enumerate(allpowers)))]
            difference = sub(Anu[N],ext(Bnu[N]))
            check(difference==rhs,'free-killed Abel absorption identity')
            direct = scale(sub(sub(Anu[N],Anu[N+1]),ext(sub(Bnu[N],Bnu[N+1]))),N)
            derivative = scale(sub([sum(items) for items in zip(*(sub(row[N-j],row[N-j+1]) for j,row in enumerate(allpowers)))],mv(A,rho[N])),N)
            check(direct==derivative,'telescoped Abel derivative')
            internal_only = [r[:d]+[F()]*(m-d) for r in rho[:N]]
            wrong = [sum(items) for items in zip(*(powers(A,internal_only[j],N-j)[-1] for j in range(N)))]
            defect = sub(difference,wrong)
            check(all(x>=0 for x in defect) and any(x>0 for x in defect),'dropping exterior flux is invalid')
            exact.append({'s':str(s),'N':N,'difference':receipt(difference),'derivative':receipt(direct),
                          'external_omission_defect_mass':str(sum(defect))})
    lf,kf = np.array(L,float),np.array(K,float)
    lam,U = np.linalg.eigh(lf); kap,Z = np.linalg.eigh(kf)
    check(abs(lam[0])<1e-12 and min(lam[1:])>0,'connected graph numerical zero eigenvalue')
    lam[0]=0.
    nf = np.array(nu,float)
    cn,ck = U.T@nf,Z.T@nf[:d]
    def integrate(N,order):
        x,w = np.polynomial.legendre.leggauss(order)
        edges = np.linspace(-24*math.log(2),24*math.log(2),17)
        vs = np.concatenate([(a+b)/2+(b-a)*x/2 for a,b in zip(edges[:-1],edges[1:])])
        weights = np.concatenate([(b-a)*w/2 for a,b in zip(edges[:-1],edges[1:])])
        ss=np.exp(vs)[:,None]
        aa=ss/(ss+lam); bb=ss/(ss+kap)
        first=(N*aa**N*(lam/(ss+lam))*cn)@U.T
        second=(N*bb**N*(kap/(ss+kap))*ck)@Z.T
        first[:,:d]-=second
        return float(weights@(first*first).sum(axis=1))
    norm2=sum(x*x for x in nu)
    gapL=m*min(J[i][j] for i in range(m) for j in range(m) if i!=j)
    gapK=min(V)
    maxL=2*max(degree); maxK=maxL+max(V)
    numeric=[]
    for N in (1,2,4,8,16):
        coarse,fine=integrate(N,qorder),integrate(N,2*qorder)
        error=abs(coarse-fine)
        check(error<=1e-8*max(1.,fine),'registered quadrature refinement')
        lower_tail=F(N*N)*norm2*(1/gapL+1/gapK)**2/F(2**49)
        upper_tail=F(N*N)*norm2*(maxL+maxK)**2/F(2**49)
        ratio=(fine+float(lower_tail+upper_tail))/(math.sqrt(N)*float(sum(nu)))
        check(ratio<280,'analytic-order safe constant numerical diagnostic')
        numeric.append({'N':N,'coarse':coarse,'fine':fine,'refinement_difference':error,
                        'small_s_tail_bound':str(lower_tail),'large_s_tail_bound':str(upper_tail),
                        'energy_plus_tail_over_sqrtN_cap_mass':ratio})
    return {'states':m,'inside':d,'whole_mass':str(sum(nu)),'exterior_absorption_mass':str(sum(mu[d:])),
            'source':receipt(nu),'absorption':receipt(mu),'exact':exact,'numeric':numeric,
            'checks':COUNT-before,'status':'PASS'}

reg=json.loads(REG.read_text())
rounds=[run(r['states'],r['inside'],r['quadrature_order_per_panel']) for r in reg['rounds']]
result={'status':'PASS','rounds':rounds,'exact_and_numeric_predicates':COUNT,
        'scope':'fixed symmetric jump operator ingredients, not actual cube/FIRST or a dimension-uniform numerical proof',
        'quadrature_error':'refinement diagnostic only, not interval enclosure; analytic omitted tails separately recorded',
        'registration_sha256':hashlib.sha256(REG.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out=BASE/'nonlocal_abel_exact_numeric_guard_results_20261007.json'
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':COUNT,'output':str(out),'max_refinement_difference':max(x['refinement_difference'] for r in rounds for x in r['numeric'])}))
