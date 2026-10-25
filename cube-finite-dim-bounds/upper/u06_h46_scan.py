"""Finite H4+H6 cost discovery. All radial real roots/positive intervals are used.
No floating sphere/cube extrema or angular integral is a cost certificate.
"""
import hashlib,json,math
from pathlib import Path
from itertools import combinations_with_replacement
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import roots_jacobi
from scipy.optimize import minimize,minimize_scalar

B4S=(0.,.002,.004,.006,.008)
B6S=(-.001,-.0005,0.,.0005,.001,.002)
class Found(Exception):
    def __init__(self,row):self.row=row

def hs(y,n):
    y=np.asarray(y);s=y.sum(axis=-1);s2=(y*y).sum(axis=-1);s3=(y*y*y).sum(axis=-1)
    return s2-3*s*s/(n+2),s3-15*s*s2/(n+8)+30*s**3/((n+8)*(n+4))

def sphere_extrema(n,B4,B6,R):
    d=B6*R**6;e=B4*R**4-15*d/(n+8);c=-3*B4*R**4/(n+2)+30*d/((n+8)*(n+4))
    vals=[]
    for m in range(1,n+1):
        vals.append((d/m**2+e/m+c,{'support':m,'uniform':True}))
        if d==0:continue
        T=-2*e/(3*d)
        for k in range(1,m):
            ell=m-k
            if k==ell:continue # degenerate critical continuum is constant, uniform represents it
            u=(1-ell*T)/(k-ell);v=T-u
            if min(u,v)>0 and max(u,v)<1:
                val=d*(k*u**3+ell*v**3)+e*(k*u*u+ell*v*v)+c
                vals.append((val,{'support':m,'k':k,'u':u,'v':v}))
    return min(vals,key=lambda z:z[0]),max(vals,key=lambda z:z[0])

@lru_cache(maxsize=2048)
def cube_min(n,B4,B6):
    mesh=8 if n==3 else 4
    tails=np.array(list(combinations_with_replacement(np.arange(mesh+1)/mesh,n-1)))[:,::-1]
    def f(z):
        y=np.r_[1.,np.asarray(z)**2];h4,h6=hs(y,n)
        return float(y.sum()**((2-n)/2)-B4*h4-B6*h6)
    values=np.array([f(z) for z in tails]);i=int(np.argmin(values));best=(float(values[i]),tails[i].tolist());calls=0
    for j in np.argsort(values)[:3]:
        opt=minimize(f,tails[j],method='L-BFGS-B',bounds=[(0.,1.)]*(n-1),options={'maxfun':80,'maxiter':30,'ftol':1e-13})
        calls+=int(opt.nfev)
        if opt.fun<best[0]:best=(float(opt.fun),opt.x.tolist())
    return best[0],{'face_point':[1.]+best[1],'grid_points':len(tails),'local_nfev':calls}

@lru_cache(maxsize=8)
def sphere_rule(n,fine=False):
    if n==3:
        nu,nphi=(24,48) if fine else (12,24)
        z,w=leggauss(nu);u=(z+1)/2;wu=w/2
        z,w=leggauss(nphi);phi=(z+1)*np.pi/4;wf=w*np.pi/4
        U,P=np.meshgrid(u,phi,indexing='ij')
        y=np.array([(1-U*U)*np.cos(P)**2,(1-U*U)*np.sin(P)**2,U*U]).reshape(n,-1).T
        weights=(8*wu[:,None]*wf[None,:]).ravel()
    else:
        nu,nv,nphi=(12,12,24) if fine else (8,8,16)
        z,w=roots_jacobi(nu,.5,0);u=(z+1)/2;wu=w/(2*np.sqrt(2))*np.sqrt(1+u)
        z,w=leggauss(nv);v=(z+1)/2;wv=w/2
        z,w=leggauss(nphi);phi=(z+1)*np.pi/4;wf=w*np.pi/4
        U,V,P=np.meshgrid(u,v,phi,indexing='ij')
        y=np.array([(1-U*U)*(1-V*V)*np.cos(P)**2,(1-U*U)*(1-V*V)*np.sin(P)**2,(1-U*U)*V*V,U*U]).reshape(n,-1).T
        weights=(16*wu[:,None,None]*wv[None,:,None]*wf[None,None,:]).ravel()
    h4,h6=hs(y,n)
    return y,h4,h6,weights

def radial_intervals(n,A,B4,B6,R,h4,h6):
    if n==3:coef=np.array([-B6*h6*R**7,0.,-B4*h4*R**5,0.,0.,0.,-A*R,1.])
    else:coef=np.array([-B6*h6*R**8,-B4*h4*R**6,0.,-A*R**2,1.])
    coef=np.trim_zeros(coef,'f');rr=np.roots(coef)
    real=sorted(float(z.real) for z in rr if abs(z.imag)<1e-7 and 0<z.real<1)
    cuts=[0.]+real+[1.];intervals=[];near=False
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        if hi-lo<1e-8:near=True
        if np.polyval(coef,(lo+hi)/2)>0:
            a,b=(R*lo,R*hi) if n==3 else (R*np.sqrt(lo),R*np.sqrt(hi))
            if b>a:intervals.append((float(a),float(b)))
    def J(r):return r*r/2-A*r**n/n-B4*h4*r**(n+4)/(n+4)-B6*h6*r**(n+6)/(n+6)
    value=sum(J(b)-J(a) for a,b in intervals)
    return value,intervals,near

def mass(n,A,B4,B6,R,fine=False):
    y,h4,h6,w=sphere_rule(n,fine);values=[];maximum=0;witness=None;near=0
    for j in range(len(w)):
        val,ints,flag=radial_intervals(n,A,B4,B6,R,float(h4[j]),float(h6[j]))
        values.append(val);near+=int(flag)
        if len(ints)>maximum:
            maximum=len(ints)
            if maximum>1:witness={'squared_direction':y[j].tolist(),'h4':float(h4[j]),'h6':float(h6[j]),'positive_intervals':ints}
    return float(w@np.array(values)),{'sphere_nodes':len(w),'maximum_positive_intervals':maximum,'multiband_witness':witness,'near_repeated_root_directions':near}

def evaluate(n,A,B4,B6,fine=False):
    R=3.2 if n==3 else 3.6
    if not(-.002<=B4<=.012 and -.003<=B6<=.003):return None
    mn,mx=sphere_extrema(n,B4,B6,R);outer_level=R**(2-n)-mn[0]
    m,location=cube_min(n,B4,B6)
    if A<.05 or A<=outer_level+1e-7 or A>=m-1e-7:return None
    M,info=mass(n,A,B4,B6,R,fine)
    row={'n':n,'A':float(A),'B4':float(B4),'B6':float(B6),'R':R,
         'outer_negative_margin':float(A-outer_level),'sphere_correction_min':float(mn[0]),'sphere_min_pattern':mn[1],
         'cube_min_before_A_numeric':m,'cube_positive_margin_numeric':m-A,'cube_min_numeric_location':location,
         'mass_dimensionless':M,'cost':M/(2**n*(m-A)),**info}
    if row['cost']<(8 if n==3 else 16):raise Found(row)
    return row

def optimize_A(n,B4,B6):
    R=3.2 if n==3 else 3.6
    mn,_=sphere_extrema(n,B4,B6,R);lower=max(.05,R**(2-n)-mn[0]);m,loc=cube_min(n,B4,B6)
    if m-lower<=2e-6:return {'n':n,'B4':B4,'B6':B6,'status':'discovery feasible interval empty','lower':lower,'cube_numeric':m}
    cache={}
    def objective(A):
        row=evaluate(n,float(A),B4,B6);cache[float(A)]=row
        return row['cost'] if row else 1e6
    opt=minimize_scalar(objective,bounds=(lower+1e-6,m-1e-6),method='bounded',options={'maxiter':40,'xatol':1e-8})
    row=cache[float(opt.x)];row.update({'optimizer_nfev':int(opt.nfev),'optimizer_success':bool(opt.success)})
    return row

def main():
    rows=[];polish=None;early=None;cross=[];refined=[];cross_refined=[]
    try:
        for b4 in B4S:
            for b6 in B6S:rows.append(optimize_A(3,b4,b6))
        best=min((z for z in rows if 'cost'in z),key=lambda z:z['cost'])
        def joint(v):
            z=evaluate(3,*map(float,v))
            return z['cost'] if z else 1e4+sum(np.square(v))
        opt=minimize(joint,[best['A'],best['B4'],best['B6']],method='Nelder-Mead',options={'maxfev':120,'xatol':1e-7,'fatol':1e-7})
        polish=evaluate(3,*map(float,opt.x));polish.update({'joint_nfev':int(opt.nfev),'joint_success':bool(opt.success)})
        ranked=sorted([z for z in rows if 'cost'in z]+[polish],key=lambda z:z['cost'])
        refined=[evaluate(3,z['A'],z['B4'],z['B6'],True) for z in ranked[:2]]
        for z in ranked[:2]:cross.append(optimize_A(4,z['B4'],z['B6']))
        valid=[z for z in cross if 'cost'in z]
        if valid:
            z=min(valid,key=lambda z:z['cost']);cross_refined=[evaluate(4,z['A'],z['B4'],z['B6'],True)]
    except Found as event:early=event.row
    data={'batch':'U06','normalization':'actual F=c_n times dimensionless F; W=1',
          'status':'floating discovery; no continuous cost certificate','rows':rows,'polish':polish,
          'n3_refinements':refined,'n4_cross':cross,'n4_refinement':cross_refined,
          'threshold_candidate':early,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    p=Path(__file__).with_name('u06_results.json');p.write_text(json.dumps(data,indent=2)+'\n')
    summary={k:data[k] for k in ('batch','polish','n3_refinements','n4_cross','n4_refinement','threshold_candidate','code_sha256')}
    summary['feasible_grid_pairs']=sum('cost'in z for z in rows)
    summary['completed_grid_pairs']=len(rows)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
