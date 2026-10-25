"""Finite exact checks for the locked U10 conditional two-basis interface.

This is not a source-domain coverage or a feasible-weight certificate.
"""
from pathlib import Path
from itertools import combinations
import argparse,hashlib,json,time
import sympy as s
from flint import arb,ctx

HERE=Path(__file__).resolve().parent

def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]

def run(root):
    started=time.monotonic()
    sources=json.loads((HERE/'029_source_snapshot.json').read_text())
    for name,h in sources.items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source drift')
    V=[s.Matrix(x) for x in [(1,-1,-1),(1,-1,1),(1,1,-1),(1,1,1)]]
    bases=[]
    for ids in combinations(range(4),3):
        M=s.Matrix.hstack(*(V[i] for i in ids));inv=M.inv()
        row_norms=[sum(inv[i,j]**2 for j in range(3)) for i in range(3)]
        remaining=next(i for i in range(4) if i not in ids)
        coords=inv*V[remaining]
        assert abs(M.det())==4 and all(v<=1 for v in row_norms)
        assert all(v in (-1,1) for v in coords)
        bases.append({'columns':ids,'determinant':str(M.det()),
                      'inverse_row_norm_squared':list(map(str,row_norms)),
                      'remaining_coordinates':list(map(str,coords))})
    centers=[sign*s.Rational(3,4)*s.eye(3)[:,i] for i in range(3) for sign in (-1,1)]
    distances=[s.expand((a-b).dot(a-b)) for a,b in combinations(centers,2)]
    assert sorted(distances)==[s.Rational(9,8)]*12+[s.Rational(9,4)]*3
    assert min(distances)>1
    a=s.Rational(5,4);buffer=s.Rational(3,2)-a
    assert a==s.Rational(3,4)+s.Rational(1,2) and buffer==s.Rational(1,4)
    # Reconstruct the radial mass and scaling, not author's code.
    r=s.symbols('r',nonnegative=True)
    mass0=s.integrate(4*s.pi*r*r*(1-r*r)**2,(r,0,1));mass1=6*s.Rational(1,2)**3*mass0
    assert mass0==32*s.pi/105 and mass1==8*s.pi/35 and mass1/mass0==s.Rational(3,4)
    # Hessian eigenvalues inside unit ball have absolute value <=8.
    z=s.symbols('z',nonnegative=True)
    radial=-4+12*z;tangential=-4+4*z
    assert [radial.subs(z,v) for v in (0,1)]==[-4,8]
    assert [tangential.subs(z,v) for v in (0,1)]==[-4,0]
    assert 8/s.Rational(1,2)**2==32
    # All nearline constants: determinant 4, event probability 1/2,
    # transverse factors (5*pi*u)^-1, longitudinal factor u/(pi*D).
    u,y=s.symbols('u y',positive=True)
    probability=s.integrate(u/(s.pi*(u*u+y*y)),(y,-u,u))
    assert probability==s.Rational(1,2)
    A,rho=s.symbols('A rho',positive=True)
    denom=1+(A+2)**2
    near_density=s.factor(s.Rational(1,4)*probability*(1/(5*s.pi*u))**2*u/(s.pi*denom))
    assert s.simplify(near_density-1/(200*s.pi**3*denom*u))==0
    near_integral=s.integrate(near_density,(u,rho,1))
    assert s.simplify(near_integral-s.log(1/rho)/(200*s.pi**3*denom))==0
    # Pair-law estimate: two restricted densities each <=t/(pi*delta^2),
    # a third unrestricted one <=1/(pi*t), all divided by det4.
    t,delta,d=s.symbols('t delta d',positive=True)
    pair_density=s.factor((t/(s.pi*delta**2))**2/(s.pi*t)/4)
    assert pair_density==t/(4*s.pi**3*delta**4)
    small=6*pair_density/t;large=1/(4*s.pi**3)
    cutoff=(large/small)**s.Rational(1,4)
    integrated=s.simplify(small*cutoff**2/2+large/(2*cutoff**2))
    constant=s.simplify(integrated.subs(delta,d/(8*s.sqrt(3)))*d*d)
    assert s.simplify(constant-48*s.sqrt(6)/s.pi**3)==0
    assert 3*s.sqrt(3)*d/(8*s.sqrt(3))<d
    # Distances to every full frame line, not just the two selected lines.
    tau=s.symbols('tau',real=True);eta=s.symbols('eta',nonnegative=True)
    v=s.Matrix([1,1,1]);w=s.Matrix([1,1,-1]);p=v+tau*w
    dist2=[s.expand(p.dot(p)-p.dot(q)**2/3) for q in V]
    for poly in dist2:
        shifted=s.Poly(s.expand((poly-s.Rational(8,3)).subs(tau,eta+1)),eta)
        assert all(co>=0 for co in shifted.all_coeffs())
    assert s.Rational(8,3)>(1+a/s.Rational(5,2))**2 # perturbation still >1
    # Exact coarse constants used for t>=60: log2>2/3, sqrt6<5/2.
    assert s.Rational(6)<s.Rational(25,4)
    upper_coarse=48*s.Rational(5,2)/60**2
    c_half_coarse=(4*s.Rational(2,3)-1)/48
    coarse_gap=c_half_coarse-upper_coarse
    assert coarse_gap==s.Rational(1,720)>0
    jump_lower=s.integrate(s.Rational(1,2)/y**2,(y,t,2*t))/s.pi
    assert jump_lower==1/(4*s.pi*t) # multiply by c
    # Four directions; line section length <=2a/sqrt3; |s|>=t/2.
    lb_constant=s.simplify(4*(2*a/s.sqrt(3))*4/s.pi)
    assert lb_constant==40/(s.sqrt(3)*s.pi)
    W,c=s.symbols('W c',nonnegative=True)
    source_gap=s.factor(c/(4*s.pi*t)-40*W/(s.sqrt(3)*s.pi*t*t))
    assert s.simplify(source_gap-(s.sqrt(3)*c*t-160*W)/(4*s.pi*s.sqrt(3)*t*t))==0
    # Strict support margins use exact elementary log2 inequalities; Arb also
    # checks their actual real values without floating comparisons.
    ctx.prec=192;pi=arb.pi();T0=4*arb(2).log()-1;cball=T0/(24*pi**3)
    jensen=arb(5)/(16*pi**3)
    margin_buffer=jensen/(arb(9)/4)-cball
    margin_support=jensen/(arb(25)/16)-cball
    tail_gap=cball/2-48*arb(6).sqrt()/(pi**3*60**2)
    assert arb(2).log()<1 and arb(2).log()>arb(2)/3
    assert margin_buffer>1/(72*pi**3)
    assert margin_support>3/(40*pi**3)
    assert tail_gap>1/(720*pi**3)
    assert arb(5)/4-(arb(5)/4)/arb(3).sqrt()>0
    return {'status':'PASS exact finite constants for U10 conditional all-space interface; no feasible weights certified',
            'source_sha256':sources,'basis_checks':bases,
            'small_ball_pair_distance_squared_counts':{'9/8':12,'9/4':3},
            'support_radius':str(a),'support_buffer_at_least':str(buffer),
            'masses':{'b0':str(mass0),'b1':str(mass1),'ratio':str(mass1/mass0)},
            'gradient_Lipschitz_bounds':{'b0':8,'b1_disjoint_six_balls':32},
            'nearline_density_lower':str(near_density),'nearline_G_lower':str(near_integral),
            'pair_density_upper':str(pair_density),'integrated_G_upper':str(constant/d**2),
            'all_four_transverse_distance_squared':list(map(str,dist2)),
            'transverse_min_distance_squared_on_tau_1_2':'8/3',
            'coarse_c_half_minus_tail_at_60_times_pi_cubed_lower':str(coarse_gap),
            'mu_lower_for_t_ge_60':'c/(4*pi*t)','LB_upper_for_t_ge_5_2':str(lb_constant*W/t**2),
            'sufficient_far_t':'max(60,160*W/(sqrt(3)*c))',
            'sufficient_far_radius':'sqrt(3)*max(60,160*W/(sqrt(3)*c))+5/4',
            'actual_strict_margins':{'G_minus_c_at_radius_3_2_interval':ends(margin_buffer),
                                     'G_minus_c_at_support_radius_interval':ends(margin_support),
                                     'c_half_minus_G_upper_at_t60_interval':ends(tail_gap)},
            'proof_only_steps':'C1,1 -> continuous L1 generator, Omega-exterior settlement using F>=0, zero-volume level set, finite measure and all-scale comparison audited in029_audit.md',
            'shared_libraries':{'SymPy':s.__version__,'Arb':'python-flint192bit','prior_proofs':'026 true frame/transfer;027 Jensen; no repeat'},
            'scope':'fixed original 3D frame and two specific bump bases; conditional sufficient interface only, no search/LP/source coverage/feasible weights/cost bound; cb10 excluded',
            'seconds':time.monotonic()-started}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args()
    result=run(args.project);result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:(HERE/'029_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
