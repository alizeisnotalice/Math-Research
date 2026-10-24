#!/usr/bin/env python3
"""Three fixed representative checks, not a repeat of the 20-box comparison.

Recompute smooth center/gradient and singular lower terms independently.
Saved full-box derivatives are used only for checking reported arithmetic.
Check eight production polynomial shifts and forty exact partial derivatives
against independently constructed SymPy coefficients. No optimizer or grid.
"""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import argparse
import importlib.util
import json
import sys
import sympy as sp

sys.dont_write_bytecode = True
ROOT=Path(__file__).resolve().parent
R=Q(15,8)
E={2:Q(69,100),6:Q(19,50),7:Q(11,10),8:Q(3,4)}


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=obj
    spec.loader.exec_module(obj)
    return obj


def verify(core_path):
    data=json.loads((ROOT/'comparison.json').read_text())
    source=ROOT/'compare_bounds.py'
    if sha256(source.read_bytes()).hexdigest()!=data['comparison_source_sha256'] or sha256(core_path.read_bytes()).hexdigest()!=data['core_sha256']:
        raise ValueError('source or reviewed core hash mismatch')
    if sha256((ROOT/'centered_bounds.md').read_bytes()).hexdigest()!=data['formula_sha256']:
        raise ValueError('formula hash mismatch')
    c=load('round017_three_center_repro_log',core_path)
    r,q=sp.symbols('r q')
    u,v,s=r*(1+q)/2,r*(1-q)/2,1-r
    a,b,zp,zm=1-v,1-u,1+r*q,1-r*q
    polynomials={}
    for p in E:
        n=p-1
        P=sp.Poly(sp.expand(sum((a**(n-k)*(zp**k+s**k)+b**(n-k)*(zm**k+s**k))/sp.Integer(n-k) for k in range(n))),r,q)
        C=sp.Poly(sp.expand(zp**n+s**n),r,q)
        polynomials[p]=(P,C)
    rs=sp.Rational(15,8)
    T=sp.Poly(sp.expand(2/(rs-1)*(2*sum(r**j/(sp.Integer(j+2)*rs**(j+1)) for j in range(29))+2*sum((r*q)**(2*k)/(sp.Integer(2*k+2)*rs**(2*k+1)) for k in range(15)))),r,q)
    def ev(poly,x,y):
        return sum((Q(int(co.p),int(co.q))*x**i*y**j for (i,j),co in poly.terms()),Q(0))
    def wlo(l,h,n):
        if h==0:
            return Q(0)
        if l==0:
            return h**n*c.log_point(h).lo if h<=Q(1,4) else -Q(1,2*n)
        return min(x**n*y for x in (l,h) for y in (c.log_point(l).lo,c.log_point(h).hi))
    def absolute(pair):
        return max(abs(Q(pair[0])),abs(Q(pair[1])))
    checks=[]
    for index in (0,12,16):
        row=data['boxes'][index]
        rl,rh=map(Q,row['r']); ql,qh=map(Q,row['q'])
        rc,qc=(rl+rh)/2,(ql+qh)/2
        rr,rq=(rh-rl)/2,(qh-ql)/2
        uu,zz=rc*(1+qc)/2,1+rc*qc
        lu,lz=c.log_point(uu),c.log_point(zz)
        du,dz=[(1+qc)/2,rc/2],[qc,rc]
        A=c.I.point(ev(T,rc,qc))
        gradient=[c.I.point(ev(T.diff(var),rc,qc)) for var in (r,q)]
        for p,(P,C) in polynomials.items():
            n,alpha=p-1,2*p*E[p]
            pp,cc=ev(P,rc,qc),ev(C,rc,qc)
            A-=alpha*(pp+cc*lu-zz**n*lz)
            for j,var in enumerate((r,q)):
                gp,gc=ev(P.diff(var),rc,qc),ev(C.diff(var),rc,qc)
                gradient[j]-=alpha*(gp+gc*lu+cc*du[j]/uu-(n*zz**(n-1)*dz[j]*lz+zz**n*dz[j]/zz))
        if A.pair()!=row['smooth_center_value'] or [g.pair() for g in gradient]!=row['smooth_center_gradient']:
            raise ValueError('independent center or gradient mismatch')
        sl,sh,zl,zh=1-rh,1-rl,1-rh*qh,1-rl*ql
        vh=rh*(1-ql)/2
        cl=sum((2*p*e*(zl**(p-1)+sl**(p-1)) for p,e in E.items()),Q(0))
        favorable=-cl*min(Q(0),c.log_point(vh).hi) if cl and vh else Q(0)
        singular=favorable+sum((2*p*e*(wlo(zl,zh,p-1)+2*wlo(sl,sh,p-1)) for p,e in E.items()),Q(0))
        if str(singular)!=row['singular_remainder']['total_lower']:
            raise ValueError('independent singular bound mismatch')
        mean=A.lo-absolute(row['smooth_uniform_gradient'][0])*rr-absolute(row['smooth_uniform_gradient'][1])*rq
        H=row['smooth_uniform_hessian']
        second=A.lo-absolute(row['smooth_center_gradient'][0])*rr-absolute(row['smooth_center_gradient'][1])*rq-(absolute(H[0][0])*rr**2+2*absolute(H[0][1])*rr*rq+absolute(H[1][1])*rq**2)/2
        if (str(mean),str(second),str(mean+singular),str(second+singular))!=(row['smooth_mean_lower'],row['smooth_second_order_lower'],row['new_mean_lower'],row['new_second_order_lower']):
            raise ValueError('saved derivative-bound arithmetic mismatch')
        checks.append({'box_index':index,'label':row['label'],
                       'center_exactly_reproduced':True,'point_gradient_exactly_reproduced':True,
                       'singular_lower_exactly_reproduced':True,
                       'saved_uniform_derivative_bounds_used_only_for_arithmetic_check':True,
                       'new_best_display':float(max(mean,second)+singular)})
    # Import only the implementation's Poly class; do not call main/prepare/build_A.
    production=load('round017_polynomial_repro_only',source)
    class Noop:
        def check(self):
            pass
    shift_checks=derivative_checks=0
    rc,qc=Q(109,250),Q(1,2000)
    def coefficients(poly):
        return {ij:Q(int(co.p),int(co.q)) for ij,co in poly.terms() if co}
    for pair in polynomials.values():
        for polynomial in pair:
            prod=production.Poly(coefficients(polynomial))
            shifted=prod.shift(rc,qc,Noop())
            wanted=sp.Poly(sp.expand(polynomial.as_expr().subs({r:r+sp.Rational(rc.numerator,rc.denominator),q:q+sp.Rational(qc.numerator,qc.denominator)},simultaneous=True)),r,q)
            if shifted.d!=coefficients(wanted):
                raise ValueError('polynomial shift mismatch')
            shift_checks+=1
            for axes,reference in [([0],polynomial.diff(r)),([1],polynomial.diff(q)),([0,0],polynomial.diff(r).diff(r)),([0,1],polynomial.diff(r).diff(q)),([1,1],polynomial.diff(q).diff(q))]:
                obj=prod
                for axis in axes:
                    obj=obj.diff(axis)
                if obj.d!=coefficients(reference):
                    raise ValueError('polynomial derivative mismatch')
                derivative_checks+=1
    # Metadata checks on all saved rows do not recompute their functions.
    if len(data['boxes'])!=20 or len({tuple(x['r']+x['q']) for x in data['boxes']})!=20:
        raise ValueError('saved box count')
    summary={'old_nonnegative_boxes':sum(Q(x['old_natural']['lower'])>=0 for x in data['boxes']),
             'new_nonnegative_boxes':sum(Q(x['new_best_lower'])>=0 for x in data['boxes']),
             'new_strictly_tighter_boxes':sum(Q(x['new_minus_old'])>0 for x in data['boxes'])}
    if summary!=data['summary']:
        raise ValueError('summary counts')
    diagnostic=json.loads((ROOT/'active_center_diagnostic.json').read_text())
    if not (diagnostic['additional_point_bound_evaluations']==1 and diagnostic['box_index']==12 and diagnostic['center']==data['boxes'][12]['center'] and Q(diagnostic['lower'])>0 and float(Q(diagnostic['lower']))==diagnostic['lower_display']):
        raise ValueError('saved extra diagnostic metadata/sign')
    return {'status':'three_representative_centers_gradients_singular_terms_and_polynomial_code_verified',
            'execution':'portable_verifier','boxes_rechecked':[0,12,16],'checks':checks,
            'exact_polynomial_shift_checks':shift_checks,'exact_polynomial_partial_derivative_checks':derivative_checks,
            'comparison_rounds_rerun':0,'uniform_derivative_enclosures_independently_recomputed':0,
            'summary_counts_verified':summary,
            'additional_center_G_diagnostic':'saved exact sign and metadata checked; not rerun',
            'comparison_source_sha256':data['comparison_source_sha256'],
            'comparison_result_sha256':sha256((ROOT/'comparison.json').read_bytes()).hexdigest(),
            'verifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Three smooth centers/gradients/singular terms and exact polynomial implementation checks only; not an independent 20-box Hessian recomputation or a global certificate.'}


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--core',type=Path)
    args=parser.parse_args()
    paths=[args.core] if args.core else [ROOT.parent/'round008/certify.py',ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments/round008/certify.py']
    core=next((p for p in paths if p.is_file()),None)
    if core is None:
        raise FileNotFoundError('Pass --core <reviewed round008/certify.py>.')
    result=verify(core)
    (ROOT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
