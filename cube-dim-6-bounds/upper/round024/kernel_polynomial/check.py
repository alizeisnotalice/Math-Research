#!/usr/bin/env python3
"""Complete exact saved N/D/P and positive-root checker; no kernel calls."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,json,math,os,time
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as QP,fmpz_poly as ZP
from sympy import Poly,Symbol,ZZ,oo

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[2]/'upper/round023/partial_fraction_kernel'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def qp(obj):
    denominator=int(obj['positive_common_denominator']);assert denominator>0
    return QP([Q(int(x),denominator)for x in obj['integer_coefficients_ascending']])

def zi(coeffs):return ZP([int(x)for x in coeffs])

def ev(poly,q):return QP([Q(int(x))for x in poly.coeffs()])(Q(q.numerator,q.denominator))

def variation(values):
    last=0;count=0
    for value in values:
        sign=(value>0)-(value<0)
        if sign:
            if last and sign!=last:count+=1
            last=sign
    return count

def dyadics(obj):
    a=F(int(obj['lower_mantissa']))*F(2)**int(obj['lower_binary_exponent'])
    b=F(int(obj['upper_mantissa']))*F(2)**int(obj['upper_binary_exponent'])
    assert a<=b
    return a,b

def saved(deadline=None):
    tic=time.perf_counter()
    def guard():
        if deadline is not None:assert time.time()<deadline,'checker time gate'
    manifest=json.loads((HERE/'input_manifest.json').read_text())
    for row in manifest['inputs']:assert sha(row['path'])==row['sha256']
    constants=json.loads((HERE/'constants.json').read_text())
    kernel=json.loads((HERE/'kernel_polynomials.json').read_text())
    critical=json.loads((HERE/'critical_polynomial.json').read_text())
    root_data=json.loads((HERE/'roots.json').read_text())
    table=json.loads((BASE/'partial_tables.json').read_text())
    old=json.loads((BASE/'results.json').read_text())
    assert sha(HERE/'kernel_polynomials.json')==critical['kernel_polynomials_sha256']
    assert sha(HERE/'critical_polynomial.json')==root_data['critical_polynomial_sha256']
    assert sha(BASE/'partial_tables.json')==kernel['partial_table_sha256']
    lo,hi=dyadics(constants['pi']);ell=F(3141592653,10**9);up=F(3141592654,10**9)
    assert ell<lo<=hi<up and constants['precision_bits']==512
    gamma=F(old['Gamma_G']);safe=(gamma/(128*hi**5),gamma/(128*lo**5))
    assert tuple(F(x)for x in constants['safe_G_rational_bounds_from_pi'])==safe
    gl,gu=dyadics(constants['G_from_old_Gamma'])
    assert gl>0 and gl<=safe[1] and gu>=safe[0]
    tail=F(constants['large_T_tail_bound']);assert tail>0 and max(safe[1],gu)/tail<=F(1,4100)
    factors={Q(r['s']):r['exponent']for r in kernel['D_factors']}
    assert factors=={Q(r['b_squared']):r['exponent']for r in table['denominator_factors']}
    assert len(factors)==359 and sum(factors.values())==417 and all(q>0 for q in factors)
    D=QP([1])
    for q,m in factors.items():D*=QP([q,1])**m
    N=qp(kernel['N']);assert N.degree()==414 and D.degree()==417 and N[0]==0
    assert N.gcd(D).degree()==0 and D.leading_coefficient()==1
    assert N.leading_coefficient()==Q(old['full_weight_checksum'])
    guard()
    # Independent local residues tie the synthesized N to every source slot.
    expected={(Q(r['b_squared']),r['order']):Q(r['coefficient'])for r in table['real_z_principal']}
    derivatives=[N]
    for _ in range(1,3):derivatives.append(derivatives[-1].derivative())
    checked_slots=0
    for index,(q,m)in enumerate(factors.items()):
        E,remainder=divmod(D,QP([q,1])**m);assert not remainder
        ed=[E]
        for _ in range(1,m):ed.append(ed[-1].derivative())
        ecoeff=[ed[k](-q)/math.factorial(k)for k in range(m)]
        ncoeff=[derivatives[k](-q)/math.factorial(k)for k in range(m)]
        assert ecoeff[0]!=0
        series=[]
        for k in range(m):
            value=(ncoeff[k]-sum((ecoeff[h]*series[k-h]for h in range(1,k+1)),Q(0)))/ecoeff[0]
            series.append(value)
            assert value==expected.get((q,m-k),Q(0))
            checked_slots+=1
        if index%32==0:guard()
    assert checked_slots==417
    P=N-64*Q(ell.numerator,ell.denominator)**6*Q(1,4100)*D
    Z=zi(critical['primitive_integer_P_ascending']);scale=Q(critical['positive_scale_P_equals_scale_times_primitive'])
    assert scale>0 and QP(Z)*scale==P and int(Z.content())==1
    assert Z[0]<0 and Z.leading_coefficient()<0 and Z.degree()==417<=512
    factors_sq=[(zi(row['integer_coefficients_ascending']),row['multiplicity'])for row in critical['squarefree_factors']]
    rebuild=ZP([int(critical['squarefree_content'])]);S=ZP([1])
    for f,m in factors_sq:
        assert m>=1 and f.gcd(f.derivative()).degree()==0
        rebuild*=f**m;S*=f
    assert rebuild==Z and S==zi(critical['squarefree_part_ascending'])
    assert S.gcd(S.derivative()).degree()==0 and S[0]!=0
    roots=root_data['roots'];assert len(roots)<=4
    V=variation(list(reversed(S.coeffs())))
    assert V==root_data['Descartes_variations']
    previous_z=None;previous_t=None
    for root in roots:
        a,b=F(root['z_lower']),F(root['z_upper'])
        assert 0<a<=b and (previous_z is None or previous_z<=a)
        previous_z=b
        if a==b:assert ev(S,a)==0
        else:assert ev(S,a)*ev(S,b)<0
        owners=[]
        for index,(f,m)in enumerate(factors_sq):
            present=(ev(f,a)==0)if a==b else(ev(f,a)*ev(f,b)<0)
            if present:owners.append((index,m))
        assert owners==[(root['squarefree_factor_index'],root['multiplicity'])]
        t0,t1=F(root['T_lower']),F(root['T_upper'])
        assert 0<t0<=t1 and t0*t0<=a<=b<=t1*t1
        assert previous_t is None or previous_t<=t0
        previous_t=t1
    guard()
    proof=root_data['completeness_proof']
    if proof=='Descartes_upper_equals_exact_disjoint_bracket_lower':
        assert V==len(roots)<=4
        count=V
    else:
        assert proof=='exact_Sturm_count_roots'
        sym=Poly.from_list([int(x)for x in reversed(S.coeffs())],Symbol('z'),domain=ZZ)
        count=int(sym.count_roots(0,oo))
        assert count==root_data['Sturm_count']
    assert count==len(roots)==root_data['positive_root_count']
    assert sum(r['multiplicity']for r in roots)==root_data['positive_roots_with_multiplicity']
    guard()
    result={'status':'PASS_EXACT_N_D_P_AND_COMPLETE_POSITIVE_ROOT_CERTIFICATE','pid':os.getpid(),
            'elapsed_seconds':time.perf_counter()-tic,'N_degree':N.degree(),'D_degree':D.degree(),
            'source_partial_slots_independently_reconstructed':checked_slots,
            'gcd_N_D_one':True,'N_zero_at_zero':True,'positive_P_scaling_verified':True,
            'pi_outer_bound_arithmetic_verified':True,'old_G_tail_arithmetic_verified':True,
            'positive_root_count':count,'root_completeness_proof':proof,
            'multiplicity_owners_verified':True,'positive_T_bracket_squares_verified':True,
            'new_kernel_calls':0,'new_full_G_calls':0,'source_hashes':[
                {'path':str(HERE/name),'sha256':sha(HERE/name)}for name in
                ('input_manifest.json','constants.json','kernel_polynomials.json','critical_polynomial.json','roots.json','check.py')]}
    temp=HERE/'checked_results.json.tmp';temp.write_text(json.dumps(result,indent=2)+'\n');os.replace(temp,HERE/'checked_results.json')
    print(json.dumps({k:v for k,v in result.items()if k!='source_hashes'}),flush=True)
    return result

if __name__=='__main__':saved()
