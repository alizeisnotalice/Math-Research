"""U04 finite geometry discovery: actual Delta box Green, rotated orbit average.

Finite boundary samples and non-interval quadrature: numerical evidence only.
"""
import hashlib,json,math
from pathlib import Path
from itertools import combinations_with_replacement
import numpy as np
from scipy.special import erfc
from numpy.polynomial.legendre import leggauss

def time_rule(nodes):
    z,w=leggauss(nodes);lo,hi=math.log(1e-4),math.log(20)
    s=lo+(z+1)*(hi-lo)/2
    return np.exp(s),w*(hi-lo)/2*np.exp(s)

def interval_heat(t,y,r):
    # t: (nt,), y: arbitrary points; output (nt,npoints).
    y=np.asarray(y); inside=np.abs(y)<r
    out=np.zeros((len(t),len(y)))
    small=t/r**2<.12
    if np.any(small):
        ts=t[small,None]; yy=y[None,:]
        val=sum(np.exp(-(yy-4*m*r)**2/(4*ts))-np.exp(-(yy-(4*m+2)*r)**2/(4*ts)) for m in range(-3,4))/np.sqrt(4*np.pi*ts)
        out[small]=val
    if np.any(~small):
        odd=np.arange(1,24,2)
        out[~small]=np.exp(-t[~small,None]*(np.pi*odd[None,:]/(2*r))**2)@np.cos(np.pi*odd[:,None]*y[None,:]/(2*r))/r
    out[:,~inside]=0
    return out

def survival(t,r):
    out=np.empty_like(t);small=t/r**2<.12
    ts=t[small]
    out[small]=1-2*sum((-1)**k*erfc((2*k+1)*r/(2*np.sqrt(ts))) for k in range(4))
    odd=np.arange(1,24,2);k=np.arange(12)
    out[~small]=np.exp(-t[~small,None]*(np.pi*odd[None,:]/(2*r))**2)@(4/np.pi*(-1.)**k/odd)
    return out

def chamber(n,mesh):
    tail=np.array(list(combinations_with_replacement(np.arange(mesh+1)/mesh,n-1)))[:,::-1]
    return np.c_[np.ones(len(tail)),tail]

def orbit_values(n,beta,theta,R,points,nodes):
    t,w=time_rule(nodes);r=np.array([beta,1/beta]+[1.]*(n-2))*R
    lowest_rate=float(np.sum((np.pi/(2*r))**2))
    tail_time=np.exp(-20*lowest_rate)/lowest_rate
    # On [0,1e-4], S=1 up to Gaussian exit tails (not zero).
    mass=float(1e-4+w@np.prod(np.array([survival(t,ri) for ri in r]),axis=0)+(4/np.pi)**n*tail_time)
    c,s=math.cos(theta),math.sin(theta)
    result=np.zeros(len(points))
    for i in range(n):
        for j in range(n):
            if i==j:continue
            rest=[k for k in range(n) if k not in (i,j)]
            for sign in (-1,1):
                coords=[c*points[:,i]+s*sign*points[:,j],-s*points[:,i]+c*sign*points[:,j]]+[points[:,k] for k in rest]
                heat=np.ones((len(t),len(points)))
                tail_coeff=np.ones(len(points))
                for coord,ri in zip(coords,r):
                    heat*=interval_heat(t,coord,ri)
                    tail_coeff*=np.where(np.abs(coord)<ri,np.cos(np.pi*coord/(2*ri))/ri,0)
                result+=(w@heat+tail_time*tail_coeff)/(2*n*(n-1))
    return mass,result

def evaluate(n,beta,theta,R,mesh,nodes):
    pts=chamber(n,mesh);mass,g=orbit_values(n,beta,theta,R,pts,nodes)
    k=int(np.argmin(g));corner=np.flatnonzero(np.all(pts==1,axis=1))[0]
    return dict(n=n,beta=beta,theta=theta,R=R,mesh=mesh,nodes=nodes,
                mass=mass,sampled_min=float(g[k]),sampled_min_at=pts[k].tolist(),
                corner_value=float(g[corner]),cost_sampled=mass/(2**n*g[k]) if g[k]>0 else None,
                boundary_points=len(pts),status='sampled optimistic cost, not continuous certificate')

def main():
    rows=[evaluate(3,b,th,r,8,160) for b in (1.,1.25,1.5) for th in (0.,np.pi/8,np.pi/4) for r in (1.7,2.1,2.5)]
    best=sorted(rows,key=lambda z:z['cost_sampled'] or math.inf)[:3]
    refinements=[evaluate(3,z['beta'],z['theta'],z['R'],16,240) for z in best]
    shapes=[(z['beta'],z['theta'],z['R']) for z in best]
    baseline=(1.,0.,2.1)
    if baseline not in shapes:shapes.append(baseline)
    cross=[evaluate(4,b,th,r,4,240) for b,th,r in shapes]
    data=dict(batch='U04',rows=rows,refinements=refinements,n4_cross=cross,
              code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output=Path(__file__).with_name('u04_results.json')
    output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(best_n3=best,refinements=refinements,n4_cross=cross,results_path=str(output),code_sha256=data['code_sha256']),indent=2))

if __name__=='__main__':main()
