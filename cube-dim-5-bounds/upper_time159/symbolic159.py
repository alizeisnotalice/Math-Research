"""New finite exact time-moment identities; no physical data or scalar kernel."""
from pathlib import Path
import hashlib,json,os,sys,time
import sympy as s

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    reg,out=map(Path,sys.argv[1:]);assert not out.exists()
    r=json.loads(reg.read_bytes())
    for row in r['sources'].values():assert sha(row['path'])==row['sha256']
    t=time.monotonic();a,b,z,d,B,V=s.symbols('a b z d B V')
    E=lambda j:s.Rational(4,j+4)*(a**(j+4)-b**(j+4))/(a**4-b**4)
    Z=(a**4-b**4)/4
    alpha=(E(3)-b*E(2))/(a-b);beta=(a*E(2)-E(3))/(a-b)
    checks=[]
    for j in (2,3):assert s.cancel(s.integrate(z**(j+3),(z,b,a))/Z-E(j))==0
    checks.append('new_z2_z3_moments_exact')
    assert s.cancel(s.integrate(z**5*(z-b)/(a-b),(z,b,a))/Z-alpha)==0
    assert s.cancel(s.integrate(z**5*(a-z)/(a-b),(z,b,a))/Z-beta)==0
    assert s.cancel(alpha+beta-E(2))==0
    checks.append('positive_barycentric_coefficients_integral_identities')
    assert s.cancel(alpha.subs(b,0)-s.Rational(4,7)*a*a)==0
    assert s.cancel(beta.subs(b,0)-s.Rational(2,21)*a*a)==0
    checks.append('infinite_time_tail_coefficients')
    positive_forms={}
    for name,x in [('alpha',alpha),('beta',beta)]:
        num,den=s.fraction(s.cancel(x.subs(a,b+d)))
        for part in (num,den):
            coeff=s.Poly(part,b,d).coeffs();assert coeff and all(c>0 for c in coeff)
        positive_forms[name]=[str(s.expand(num)),str(s.expand(den))]
    checks.append('coefficientwise_positive_for_b_nonnegative_d_positive')
    theta=(E(1)-b)/(a-b)
    difference=theta*a*a*(B+V)+(1-theta)*b*b*B-alpha*(B+V)-beta*B
    num,den=s.fraction(s.cancel(difference.subs(a,b+d)))
    assert all(c>0 for c in s.Poly(num,b,d,B,V).coeffs())
    assert all(c>0 for c in s.Poly(den,b,d).coeffs())
    checks.append('old_minus_new_Taylor_coefficientwise_nonnegative')
    result={'status':'PASS_new159_exact_symbolic_identities','pid':os.getpid(),
        'seconds':time.monotonic()-t,'checks':checks,'positive_coefficient_forms':positive_forms,
        'difference_numerator':str(s.expand(num)),'difference_denominator':str(s.expand(den)),
        'real_H_R_J_factor_mean_gain_or_C5_evaluations':0,'old_baseline_replays':0,
        'registration_sha256':sha(reg)}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','pid','seconds','checks')}))
if __name__=='__main__':main()
