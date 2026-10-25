"""U09: one finite path generator + exact layer/radial mass algebra.
Not a continuum kernel/cost approximation; no new parameter search.
"""
from pathlib import Path
import hashlib,json
import sympy as s

N=5
Q=s.zeros(N)
for i in range(N-1):Q[i,i+1]=Q[i+1,i]=1
for i in range(N):Q[i,i]=-sum(Q[i,j] for j in range(N) if j!=i)
one=s.ones(N,1)
assert Q==Q.T and Q*one==s.zeros(N,1)
t=s.Rational(1,4);P=s.eye(N)+t*Q
assert all(z>=0 for z in P) and P*one==one
G=s.Matrix([0,s.Rational(3,4),s.Rational(1,2),s.Rational(1,4),0])
e=s.Matrix([0,1,0,0,0]);exit=s.Matrix([s.Rational(3,4),0,0,0,s.Rational(1,4)])
assert Q*G==exit-e
c=s.Rational(3,8)
H=G.applyfunc(lambda z:max(z-c,0))
mu=Q*H+e
assert H==s.Matrix([0,s.Rational(3,8),s.Rational(1,8),0,0])
assert mu==s.Matrix([s.Rational(3,8),s.Rational(3,8),s.Rational(1,8),s.Rational(1,8),0])
assert all(z>=0 for z in mu) and sum(mu)==1
jensen=P*H-(P*G).applyfunc(lambda z:max(z-c,0))
assert all(z>=0 for z in jensen)
assert all(z>=0 for z in P*H-H+t*e)
# The finite model has source-node mass; it is not proof of continuum mu({0})=0.

cc,gg,aa,JJ=s.symbols('c g a J',positive=True)
r,lam=s.symbols('r lam',positive=True)
rows=[]
for n in (3,4):
    alpha=n-1;q=s.Rational(n,n-1)
    rho=(aa/cc)**s.Rational(1,alpha)
    radial=s.integrate((aa*r**(-alpha)-cc)*r**(n-1),(r,0,rho))
    target=s.Rational(alpha,n)*cc**(-s.Rational(1,alpha))*aa**q
    assert s.simplify(radial-target)==0
    layer=s.integrate(lam**(-q),(lam,cc,s.oo))
    assert s.simplify(layer-alpha*cc**(-s.Rational(1,alpha)))==0
    cost=s.Rational(alpha,n)*JJ*cc**(-s.Rational(1,alpha))/(2**n*(gg-cc))
    logder=s.factor(s.diff(cost,cc)/cost)
    assert s.simplify(logder-(n*cc-gg)/(alpha*cc*(gg-cc)))==0
    opt=s.simplify(cost.subs(cc,gg/n))
    opttarget=n**s.Rational(1,alpha)*JJ/(2**n*gg**q)
    assert s.simplify(opt-opttarget)==0
    along=s.integrate(r**(2-n),(r,0,1))
    off=s.integrate(r**(4-n),(r,0,1))
    assert along==s.oo and off==s.Rational(1,5-n)
    rows.append({'n':n,'q':str(q),'radial_mass':str(target),'layer_integral':str(layer),
                 'cost_log_derivative':str(logder),'optimal_c':'g/'+str(n),
                 'optimal_mass_ratio':str(opttarget),'line_small_time_integral':str(along),
                 'off_line_small_time_integral':str(off)})
vec=lambda v:[str(x) for x in v]
data={'batch':'U09','status':'finite exact diagnostics PASS; not continuum source/cost certification',
      'matrix':{'generator':[[str(Q[i,j]) for j in range(N)] for i in range(N)],'t':str(t),
                'green':vec(G),'green_source':vec(Q*G),'threshold':str(c),'positive_part':vec(H),
                'new_positive_source':vec(mu),'new_source_mass':str(sum(mu)),
                'jensen_slack':vec(jensen),'note':'source-node positive mass3/8 exists in finite model; continuum no-atom uses small-cutoff proof'},
      'symbolic_mass_checks':rows,'new_upper_bound':False,
      'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('u09_results.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(data,indent=2,ensure_ascii=False))
