"""U05 finite cost discovery for actual F=c_n(r^(2-n)-A-B H4), W=1.

Exact edge reduction supplied by main agent; floating polynomial roots and
angular quadrature are not continuous cost certificates.
"""
import hashlib,json,math
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
from scipy.special import roots_jacobi
from scipy.optimize import minimize_scalar
import sympy as sp

BS=(0.,.001,.003,.006,.01,.016,.025,.04,.05)

def cube_min(n,b):
    candidates=[]
    t=Poly([0.,1.])
    for k in range(1,n):
        roots=[]
        if b:
            polynomial=(4*b*b/25)*(6*k-4*t)**2*(k+t)**3-1 if n==3 else b*(k-t)*(k+t)**2-1
            for z in polynomial.roots():
                if abs(z.imag)<1e-8 and 0<z.real<1:
                    x=float(z.real)
                    derivative=b*(6*k-4*x)/5-1/(2*(k+x)**1.5) if n==3 else b*(k-x)-1/(k+x)**2
                    if abs(derivative)<1e-7:roots.append(x)
        for z in [0.,1.]+roots:
            value=(k+z)**((2-n)/2)-b*(k+z*z-3*(k+z)**2/(n+2))
            candidates.append((float(value),k,float(z)))
    value,k,z=min(candidates)
    return value,{'k_ones':k,'last_square':z,'zeros':n-k-1,'value_before_A':value,'candidates':len(candidates)}

def sphere_rule(n,refine=False):
    if n==3:
        nu,nphi=(48,96) if refine else (24,48)
        zu,wu=leggauss(nu);u=(zu+1)/2;wu=wu/2
        zf,wf=leggauss(nphi);phi=(zf+1)*np.pi/4;wf=wf*np.pi/4
        U,P=np.meshgrid(u,phi,indexing='ij')
        coords=[np.sqrt(1-U*U)*np.cos(P),np.sqrt(1-U*U)*np.sin(P),U]
        weight=(8*wu[:,None]*wf[None,:]).ravel()
    else:
        nu,nv,nphi=(16,16,32) if refine else (12,12,24)
        zu,wu=roots_jacobi(nu,.5,0);u=(zu+1)/2
        wu=wu/(2*np.sqrt(2))*np.sqrt(1+u)
        zv,wv=leggauss(nv);v=(zv+1)/2;wv=wv/2
        zf,wf=leggauss(nphi);phi=(zf+1)*np.pi/4;wf=wf*np.pi/4
        U,V,P=np.meshgrid(u,v,phi,indexing='ij')
        coords=[np.sqrt(1-U*U)*np.sqrt(1-V*V)*np.cos(P),np.sqrt(1-U*U)*np.sqrt(1-V*V)*np.sin(P),np.sqrt(1-U*U)*V,U]
        weight=(16*wu[:,None,None]*wv[None,:,None]*wf[None,None,:]).ravel()
    h=(sum(c**4 for c in coords)-3/(n+2)).ravel()
    surface=4*np.pi if n==3 else 2*np.pi*np.pi
    assert abs(weight.sum()/surface-1)<1e-12
    return h,weight

def outer(n,b,A=None):
    if b==0:return 2*A**(-1/(n-2)) if A else math.inf,0.
    q=2*(n-1)/(n*(n+2))
    r=((n-2)/(4*q*b))**(1/(n+2))
    return r,(n+2)/4*r**(2-n)

def mass(n,A,b,h,w,R):
    if b==0:
        rho=np.full_like(h,A**(-1/(n-2)))
    else:
        lo=np.zeros_like(h);hi=np.full_like(h,R)
        for _ in range(60):
            mid=(lo+hi)/2
            positive=mid**(2-n)-A-b*h*mid**4>0
            lo=np.where(positive,mid,lo);hi=np.where(positive,hi,mid)
        rho=(lo+hi)/2
    radial=rho*rho/2-A*rho**n/n-b*h*rho**(n+4)/(n+4)
    assert radial.min()>0
    return float(w@radial)

def optimize(n,b):
    m,loc=cube_min(n,b);R,amin=outer(n,b)
    if (b and R<=math.sqrt(n)) or amin>=m:
        return {'n':n,'B':b,'status':'outer-shell/cube-positive conditions infeasible','R':R,'Amin':amin,'cube_min':m}
    h,w=sphere_rule(n)
    gap=m-amin
    opt=minimize_scalar(lambda A:mass(n,A,b,h,w,outer(n,b,A)[0])/(2**n*(m-A)),bounds=(amin+1e-8*gap,m-1e-8*gap),method='bounded',options={'xatol':1e-12,'maxiter':80})
    A=float(opt.x);R,amin=outer(n,b,A);M=mass(n,A,b,h,w,R)
    boundary_level=R**(2-n)+2*(n-1)/(n*(n+2))*b*R**4
    return {'n':n,'B':b,'A':A,'R':R,'Amin':amin,'outer_negative_margin':A-boundary_level,
            'cube_min_before_A':m,'cube_positive_margin':m-A,'minimum_location':loc,
            'mass_dimensionless':M,'cost':M/(2**n*(m-A)),
            'sphere_nodes':len(h),'optimizer_nfev':int(opt.nfev),'optimizer_success':bool(opt.success),
            'status':'numerical cost only; continuous source/edge reduction separate'}

def refine(row):
    n=row['n'];h,w=sphere_rule(n,True)
    M=mass(n,row['A'],row['B'],h,w,row['R'])
    return {'n':n,'B':row['B'],'A':row['A'],'sphere_nodes':len(h),
            'mass_dimensionless':M,'cost':M/(2**n*row['cube_positive_margin'])}

def main():
    harmonic=[]
    for n in (3,4):
        x=sp.symbols('x:'+str(n));r2=sum(y*y for y in x)
        H=sum(y**4 for y in x)-sp.Rational(3,n+2)*r2*r2
        lap=sp.simplify(sum(sp.diff(H,y,2) for y in x));assert lap==0
        harmonic.append({'n':n,'laplacian_H4':str(lap)})
    rows=[];early=False
    for b in BS:
        z=optimize(3,b);rows.append(z)
        if z.get('cost',math.inf)<8:early=True;break
    best=sorted([z for z in rows if 'cost'in z],key=lambda z:z['cost'])
    refined=[refine(z) for z in best[:2]]
    cross=[]
    if not early:
        for z in best[:3]:
            q=optimize(4,z['B']);cross.append(q)
            if q.get('cost',math.inf)<16:early=True;break
    valid=[z for z in cross if 'cost'in z]
    cross_refined=[refine(min(valid,key=lambda z:z['cost']))] if valid else []
    data={'batch':'U05','normalization':'F=c_n[r^(2-n)-A-BH4], W=1; c_n cancels in cost',
          'harmonic_checks':harmonic,'n3_rows':rows,'n3_refinements':refined,
          'n4_cross':cross,'n4_refinement':cross_refined,'threshold_triggered':early,
          'status':'floating discovery; no new continuous upper bound',
          'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    result=Path(__file__).with_name('u05_results.json');result.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()
