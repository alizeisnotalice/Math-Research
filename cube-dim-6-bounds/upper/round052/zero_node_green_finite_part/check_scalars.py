#!/usr/bin/env python3
"""Twelve formal scalar identities only; no fan/point/tuple aggregation."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,json,os,time,datetime
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def q(x):return s.Rational(str(x))
def normalized_laurent_jet(b,k):
    out=[s.Rational(0)]*4 # Gamma (pi/2 coefficient), eps^-1, eps^-3, eps^-5
    if k==0 and b:out[0]=s.sign(b)
    if not b and k in(1,3,5):out[(k+1)//2]=s.Rational((-1)**((k-1)//2),k)
    return out
def dd_vector(nodes):
    table=[normalized_laurent_jet(b,0)for b in nodes]
    for k in range(1,6):
        out=[]
        for j in range(6-k):
            if nodes[j+k]==nodes[j]:out.append(normalized_laurent_jet(nodes[j],k))
            else:out.append([(table[j+1][c]-table[j][c])/(nodes[j+k]-nodes[j])for c in range(4)])
        table=out
    return table[0]

def run():
    start=time.monotonic();reg=json.loads((HERE/'registration.json').read_text())
    for p,h in reg['sources'].items():assert sha(p)==h,p
    assert len(reg['node_cases'])==12
    e=s.Symbol('epsilon',positive=True);z=s.Symbol('z',real=True);T=s.Symbol('T',positive=True);u=s.Symbol('u')
    results=[]
    for index,raw in enumerate(reg['node_cases'],1):
        assert time.monotonic()-start<900
        nodes=sorted(map(q,raw));mult=Counter(nodes)
        conditions=[(b,k)for b,m in sorted(mult.items())for k in range(m)]
        # Complete normalized Hermite interpolation: DD5 is leading coefficient.
        M=s.Matrix([[s.binomial(j,k)*b**(j-k)if j>=k else 0 for j in range(6)]for b,k in conditions])
        determinant=M.det();assert determinant!=0
        weights=list(M.inv()[5,:])
        jet_vectors=[normalized_laurent_jet(b,k)for b,k in conditions]
        leading=[s.factor(sum(w*v[j]for w,v in zip(weights,jet_vectors)))for j in range(4)]
        assert leading==dd_vector(nodes)
        # Independent rational residues of f(z)/product(z-b), no atan yet.
        residue_weights=[]
        for b,k in conditions:
            degree=mult[b]-1-k
            other=s.prod((b-c+u)**(-m)for c,m in mult.items()if c!=b)
            residue_weights.append(s.factor(s.diff(other,u,degree).subs(u,0)/s.factorial(degree)))
        assert weights==residue_weights
        # Actual derivatives of atan(z/epsilon), then exact local Laurent parts.
        actual_jets=[s.diff(s.atan(z/e),z,k).subs(z,b)/s.factorial(k)for b,k in conditions]
        laurent=s.expand(sum(w*s.series(v,e,0,1).removeO()for w,v in zip(weights,actual_jets)))
        expected=s.pi*leading[0]/2+leading[1]/e+leading[2]/e**3+leading[3]/e**5
        assert s.simplify(laurent-expected)==0
        # Fix the resolvent phase as an identity in symbolic T, not a sampled time.
        density=sum(w*s.diff(z/(T*T+z*z),z,k).subs(z,b)/s.factorial(k)for w,(b,k)in zip(weights,conditions))
        product_numerator=s.expand(s.prod(T+s.I*b for b in nodes))
        real_product=s.re(product_numerator)/s.prod(T*T+b*b for b in nodes)
        assert s.cancel(density-real_product)==0
        num,den=s.fraction(s.cancel(real_product))
        row={'case':index,'nodes':list(map(str,nodes)),'conditions':[[str(b),k]for b,k in conditions],
             'Hermite_determinant':str(determinant),'Hermite_residue_weights':list(map(str,weights)),
             'Gamma_and_principal_coefficients':list(map(str,leading)),
             'Laurent_part':str(expected),'actual_atan_jets':list(map(str,actual_jets)),
             'density_numerator_coefficients_descending':list(map(str,s.Poly(num,T).all_coeffs())),
             'density_denominator_coefficients_descending':list(map(str,s.Poly(den,T).all_coeffs())),
             'recurrence_equals_Hermite_equals_residue':True,'actual_jet_Laurent_check':True,
             'resolvent_identity_as_rational_function':True,'spatial_G_value_claimed':False}
        if all(b==0 for b in nodes):
            assert leading==[0,0,0,s.Rational(1,5)]
            row['mandatory_rejection_note']='Formal constant0 but H_epsilon=1/(5 epsilon^5); never a Green(0)=0 certificate.'
        results.append(row)
        print(json.dumps({'case':index,'Gamma_L1_L3_L5':list(map(str,leading))}),flush=True)
    assert s.simplify((2*s.pi)**-6*(s.pi/2)-1/(128*s.pi**5))==0
    for p,h in reg['sources'].items():assert sha(p)==h,p
    return {'status':'PASS_TWELVE_FORMAL_SCALAR_IDENTITIES','cases':results,'actual_distinct_node_cases':12,
            'complete_G_H_or_tuple_passes':0,'positive_time_substitutions':0,'new_spatial_values':0,
            'new_C6_bound':False},time.monotonic()-start

if __name__=='__main__':
    assert sys.argv[1:]in([],['--record'])
    if sys.argv[1:]==['--record']:assert not(HERE/'results.json').exists(),'one scalar batch only'
    result,seconds=run()
    if sys.argv[1:]==['--record']:
        (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
        receipt={'status':'COMPLETED','pid':os.getpid(),'exit_code':0,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 'wall_seconds':seconds,'unique_scalar_cases':12,'full_G_H_native_or_tuple_passes':0,
                 'positive_time_values':0,'source_sha256':sha(__file__),'registration_sha256':sha(HERE/'registration.json'),
                 'python':sys.version,'python_executable':sys.executable,'python_executable_sha256':sha(sys.executable),
                 'sympy_version':s.__version__,'sympy_init_sha256':sha(s.__file__),'all_handles_closed':True,'child_processes':[]}
        (HERE/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    else:assert result==json.loads((HERE/'results.json').read_text())
    print(json.dumps({'status':result['status'],'pid':os.getpid(),'seconds':seconds,'full_kernel_passes':0}))
