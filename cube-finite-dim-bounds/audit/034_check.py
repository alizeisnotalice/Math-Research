"""Finite independent checks of whole-face short-time and convex-chord interfaces.

No Green point/node, live part, author checker, or cache calls. Full I is pending.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb
import argparse,hashlib,json,time
import sympy as s
from flint import arb,arb_mat,ctx

HERE=Path(__file__).resolve().parent
def aq(x):
    x=F(x);return arb(x.numerator)/x.denominator
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def pp(x,a,b):return x.root(b)**a
def inside(x,interval):
    lo,hi=map(F,ends(x));a,b=map(F,interval);return a<lo<=hi<b
def exact_packed(v):
    assert int(v[2])==0
    return F(int(v[0]))*F(2)**int(v[1])
def budget(n,cn,A,eps,target):
    C=aq(cn)/arb.pi()**n
    E=pp(arb(2),-(n-1)**2,n)*aq(F(n,n+1))*pp(C,1,n)*pp(aq(A)/arb.pi(),n-1,n)*pp(aq(eps),n+1,n)
    room=pp(aq(target),n-1,n)-E;assert room>0
    return E,pp(room,n,n-1)

def run(root):
    started=time.monotonic();snap=json.loads((HERE/'034_source_snapshot.json').read_text())
    for name,h in snap.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,(name,'source drift')
    primary=json.loads((root/'results/round011_face_time_budget.json').read_text())
    ctx.prec=256;models=[]
    for n,tau in [(3,F(1,3)),(4,F(0))]:
        body=[(1,)+tuple(2*((k>>j)&1)-1 for j in range(n-1)) for k in range(2**(n-1))]
        rates={v:F(1) for v in body}
        if tau:
            for i,j in combinations(range(n),2):
                for sign in (-1,1):
                    v=[0]*n;v[i]=1;v[j]=sign;rates[tuple(v)]=tau
        active=[v for v in rates if v[0]];k=len(active);A=sum((rates[v] for v in active),F(0));S=F(0);dets=[]
        for v in active:
            best=max(abs(int(s.Matrix.hstack(s.Matrix(v),*(s.Matrix(w) for w in rest)).det())) for rest in combinations(body,n-1))
            assert best=={3:4,4:16}[n];dets.append(best);S+=rates[v]/best
        cn=k*k*S;assert k==8 and cn=={3:F(256,3),4:F(32)}[n]
        q=F(n,n-1)
        assert (F(2-n)*(q-1)+1)/q==F(1,n)
        assert F(n-1)/q==F((n-1)**2,n)
        pcase=next(v for v in primary['cases'] if v['n']==n)
        assert pcase['k']==k and F(pcase['A'])==A and F(pcase['C_numerator_over_pi_n'])==cn
        records=[]
        for row in pcase['budgets']:
            E,B=budget(n,cn,A,F(row['epsilon']),F(pcase['I_target']))
            assert inside(E,row['E_interval']) and inside(B,row['allowed_I_H_interval'])
            records.append({'epsilon':row['epsilon'],'E_interval':ends(E),'allowed_IH_interval':ends(B),'strictly_inside_primary_192bit':True})
        models.append({'n':n,'tau':str(tau),'active_count':k,'active_rate_sum':str(A),'active_max_basis_dets':dets,'Cf_numerator_over_pi_n':str(cn),'time_norm_exponent':'1/'+str(n),'orthant_norm_factor_exponent':str(-F((n-1)**2,n)),'primary_budget_checks':records})
    Etime,allowed=budget(3,F(256,3),F(16,3),F(9,5000),F(29,50000))
    # The33 degree32 basis means and complete polynomial conversion identities.
    ctx.prec=384;degree=32;N=33;u=s.symbols('u')
    polys=[s.Poly(s.chebyshevt(k,2*u-1),u) for k in range(N)]
    pc=[[F(polys[k].nth(l)) for l in range(k+1)] for k in range(N)]
    bf=[[sum((pc[k][l]*F(comb(j,l),comb(degree,l)) for l in range(min(j,k)+1)),F(0)) for k in range(N)] for j in range(N)]
    means=[F(1)]+[F(0) if k%2 else F(-1,k*k-1) for k in range(1,N)]
    for k in range(N):
        assert sum((pc[k][l]/(l+1) for l in range(k+1)),F(0))==means[k]
        assert sum((bf[j][k] for j in range(N)),F(0))/N==means[k]
        back=sum((s.Rational(bf[j][k].numerator,bf[j][k].denominator)*comb(degree,j)*u**j*(1-u)**(degree-j) for j in range(N)),s.Integer(0))
        assert s.Poly(s.expand(back)-polys[k].as_expr(),u).is_zero
    A=arb_mat([[(arb.pi()*k*(2*j+1)/(2*N)).cos()*aq(F(1 if k==0 else 2,N)) for j in range(N)] for k in range(N)])
    Q=arb_mat([[(arb.pi()*k*(2*j+1)/(2*N)).cos() for k in range(N)] for j in range(N)])
    identity=A*Q
    for i in range(N):
        for j in range(N):assert identity[i,j].contains(int(i==j)) and identity[i,j].rad()<arb(2)**-300
    h=F(1,256);rho=F(5,2);eps=F(9,5000);height=h*(rho-1/rho)/2;delta=eps-height/3
    assert height==F(21,5120) and delta==F(277,640000)>0
    M=aq(F(15,1088)/delta**2)/arb.pi()**3
    error=176*M*aq(F(2,5))**32
    D=aq(F(15,1088)/((eps-eps/3)**2*eps))/arb.pi()**3
    node_error=2*D*aq(h/F(2**128))
    assert M<aq(F('2373.630728')) and error<aq(F('7.706293505929e-8')) and node_error<aq(F('3.938477182864e-36'))
    # Algebraic stable slope identity. a=sqrt(U),b=sqrt(L) positive/nonnegative.
    a,b=s.symbols('a b',positive=True)
    assert s.cancel((a*a+a*b+b*b)/(a+b)-(a**3-b**3)/(a*a-b*b))==0
    grid=128;patches=[(i,j) for i in range(grid) for j in range(i+1)]
    assert len(patches)==8256 and len(set(patches))==8256
    total_weight=sum(1 if i==j else 2 for i,j in patches);assert total_weight==grid*grid
    # All12 part index sets partition the8256 prescribed cells exactly.
    assigned=[{(i,j) for i,j in patches if (i*(i+1)//2+j)%12==p} for p in range(12)]
    assert sum(map(len,assigned))==8256 and set.union(*assigned)==set(patches)
    ctx.prec=512;pilot=json.loads((root/'upper/continuous_bridge/cb11_integral_cheb_four_result.json').read_text());chords=[]
    assert pilot['config']['epsilon']=='9/5000' and pilot['config']['degree']==32 and pilot['config']['grid']==128
    for row in pilot['patches']:
        L,U,m=map(lambda k:exact_packed(row[k]),['L','U','mean_H_upper'])
        assert 0<=L<=m<=U and U>0
        i,j=row['patch'];weight=1 if i==j else 2;assert row['weight']==weight
        slope=(aq(U)+(aq(U)*aq(L)).sqrt()+aq(L))/(aq(U).sqrt()+aq(L).sqrt())
        chord=aq(L)*aq(L).sqrt()+slope*aq(m-L)
        contribution=chord*aq(F(weight,grid*grid))
        recorded=F(row['IH_contribution_exact'])
        assert recorded==exact_packed(row['IH_contribution_upper'])
        assert F(ends(contribution)[1])<=recorded
        chords.append({'patch':[i,j],'recorded_exact_upper':str(recorded),'independent_contribution_interval':ends(contribution),'weight':weight,'scalar_bound_valid':True})
    # The full-I target remains conditional; g is independently available from033.
    gl=F(131,25000);it=F(29,50000);slack=64*gl**3-27*it**2;assert slack>0
    return {'status':'PASS whole-face short-time Lq lemma and11O conditional mean/chord/area interface; full I not independently accepted here',
            'source_sha256':snap,'short_time_models':models,
            'active_11O_budget':{'epsilon':'9/5000','I_target':'29/50000','Etime_interval':ends(Etime),'Etime_ball':Etime.str(50),'allowed_IH_interval':ends(allowed),'allowed_IH_ball':allowed.str(50)},
            'degree32_checks':{'exact_polynomial_identities':33,'exact_Cheb_and_Bernstein_mean_identities':33,'DCT_inverse_entries':1089,'ellipse_height':str(height),'ellipse_delta':str(delta),'M_interval':ends(M),'Epoly_interval':ends(error),'node_displacement_interval':ends(node_error)},
            'triangle':{'cells':8256,'total_weight':total_weight,'weighted_area':'1','all12part_partition_verified':True},
            'four_pilot_scalar_chord_checks':chords,'new_Green_or_node_evaluations':0,
            'conditional_final_target':{'g_certified_by033':'131/25000','I_not_yet_certified_target':'29/50000','exact_slack_64g3_minus27I2':str(slack)},
            'scope':'proof and finite constants/interface only; no live part files/caches read/written, no8256 nodes replay/aggregate acceptance, no new Cn bound',
            'shared_libraries':{'SymPy':s.__version__,'Fraction':'stdlib','Arb':'python-flint256/384/512bit finite constants, no full-batch precision change'},
            'seconds':time.monotonic()-started}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,default=HERE.parent);ap.add_argument('--save',action='store_true');args=ap.parse_args()
    out=run(args.project);out['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.save:(HERE/'034_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'},indent=2))
