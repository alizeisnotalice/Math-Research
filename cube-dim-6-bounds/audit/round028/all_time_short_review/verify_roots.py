#!/usr/bin/env python3
"""Independent saved rational identity/root audit; no kernel or root search."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,os,signal,sys,time
sys.set_int_max_str_digits(0)
from flint import fmpq as Q,fmpq_poly as QP,fmpz_poly as ZP,arb,ctx

HERE=Path(__file__).resolve().parent
TOP=Path('/Users/zhengzhihao/Desktop/T')
BROOT=TOP/'Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper'
B24=BROOT/'round024/kernel_polynomial'
B23=BROOT/'round023/partial_fraction_kernel'
F23=HERE.parents[1]/'round023/partial_fraction_review/results.json'
INPUTS={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):INPUTS[str(p)]=sha(p);return json.loads(Path(p).read_text())
def pin(p,h):INPUTS[str(p)]=sha(p);assert INPUTS[str(p)]==h,str(p)
def write(n,d):(HERE/n).write_text(json.dumps(d,indent=2)+'\n')
def endpoints(x):
    def cv(y):
        r=y.fmpq();return F(int(r.numerator),int(r.denominator))
    return cv(x.lower()),cv(x.upper())
def evalq(poly,x):return poly(Q(x))

def main():
    assert not (HERE/'audit_receipt.json').exists(),'one audit execution; no automatic replay'
    start=time.time();locked=read(HERE/'locked_inputs.json')
    for name,h in locked['review_sha256'].items():pin(HERE/name,h)
    for path,h in locked['frozen_sha256'].items():pin(Path(path),h)
    b28=Path(locked['B028_directory'])
    assert (b28/'report.md').exists() and (b28/'manifest.json').exists()
    for source in (B24,b28):
        manifest=read(source/'manifest.json')
        if 'files' in manifest:
            for row in manifest['files']:pin(source/row.get('path',row.get('name')),row['sha256'])
        else:
            for name,h in manifest['whitelist_sha256'].items():pin(source/name,h)
        for path,h in manifest.get('external_readonly_sha256',{}).items():pin(Path(path),h)
    receipt={'status':'INDEPENDENT_AUDIT_STARTED','pid':os.getpid(),'start_unix':start,
             'deadline_unix':start+1200,'native_H_attempts':0,'root_isolation_calls':0,'Sturm_calls':0,
             'new_PF_full_sums':0,'new_G_full_calls':0}
    write('audit_receipt.json',receipt)
    signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(TimeoutError('root audit 180s cap')))
    signal.alarm(180)
    try:
        tables=read(B23/'partial_tables.json');f23=read(F23)
        pin(F23,'a7135d22c023d1dab9251b24f282810a3e67c785e0aa6709326db7ef2933deb9')
        kernel=read(B24/'kernel_polynomials.json');critical=read(B24/'critical_polynomial.json')
        roots=read(B24/'roots.json');constants=read(B24/'constants.json')
        assert kernel['partial_table_sha256']==INPUTS[str(B23/'partial_tables.json')]
        assert critical['kernel_polynomials_sha256']==INPUTS[str(B24/'kernel_polynomials.json')]
        assert roots['critical_polynomial_sha256']==INPUTS[str(B24/'critical_polynomial.json')]
        nd=int(kernel['N']['positive_common_denominator']);assert nd>0
        N=QP([Q(int(x),nd) for x in kernel['N']['integer_coefficients_ascending']])
        Z=ZP([int(x) for x in critical['primitive_integer_P_ascending']]);Zq=QP(Z)
        scale=Q(critical['positive_scale_P_equals_scale_times_primitive']);assert scale>0
        ell,u,c=Q(3141592653,10**9),Q(3141592654,10**9),Q(1,4100)
        assert Q(critical['ell'])==ell and Q(critical['c'])==c
        factors={Q(r['s']):r['exponent'] for r in kernel['D_factors']}
        expected={(Q(r['b_squared']),r['order']):Q(r['coefficient']) for r in tables['real_z_principal']}
        assert factors=={Q(r['b_squared']):r['exponent'] for r in tables['denominator_factors']}
        assert len(factors)==len(kernel['D_factors'])==359 and sum(factors.values())==417
        assert len(expected)==len(tables['real_z_principal'])==417
        assert all(q>0 and 1<=m<=3 for q,m in factors.items())
        assert N.degree()==414 and N[0]==0 and Z.degree()==417
        assert Z[0]<0 and Z.leading_coefficient()<0 and int(Z.content())==1
        assert -scale*Z.leading_coefficient()==64*ell**6*c
        assert N.leading_coefficient()==Q(f23['full_weight_checksum'])
        nders=[N,N.derivative(),N.derivative().derivative()]
        zders=[Zq,Zq.derivative(),Zq.derivative().derivative()]
        checked=0
        for q,m in factors.items():
            # E(z)=D(z)/(z+q)^m, expanded locally by multiplying tiny
            # factor jets. No dense D is constructed or divided here.
            E=[Q(1)]+[Q(0)]*(m-1)
            for r,n in factors.items():
                if r==q:continue
                d=r-q
                f=[Q(math.comb(n,k))*d**(n-k) if k<=n else Q(0) for k in range(m)]
                E=[sum((E[h]*f[k-h] for h in range(k+1)),Q(0)) for k in range(m)]
            assert E[0]!=0
            H=[]
            for k in range(m):
                nv=nders[k](-q)/math.factorial(k)
                value=(nv-sum((E[h]*H[k-h] for h in range(1,k+1)),Q(0)))/E[0]
                H.append(value);assert value==expected.get((q,m-k),Q(0));checked+=1
                # N-scale*Z has every D zero with its exact multiplicity.
                assert nders[k](-q)==scale*zders[k](-q)
            assert H[0]!=0  # N and D have no common pole factor.
        assert checked==417
        # The 417 zeros plus the matching leading coefficient establish
        # N-scale*Z=64*ell^6*c*D exactly, without a dense D reconstruction.
        sf=[];rebuild=ZP([int(critical['squarefree_content'])]);S=ZP([1])
        for row in critical['squarefree_factors']:
            f=ZP([int(x) for x in row['integer_coefficients_ascending']]);m=row['multiplicity']
            assert m>=1 and f.gcd(f.derivative()).degree()==0
            rebuild*=f**m;S*=f;sf.append((f,m))
        assert rebuild==Z and S==ZP([int(x) for x in critical['squarefree_part_ascending']])
        assert S.gcd(S.derivative()).degree()==0 and S[0]!=0
        signs=[1 if x>0 else -1 for x in reversed(S.coeffs()) if x]
        variations=sum(a!=b for a,b in zip(signs,signs[1:]));assert variations==2
        assert roots['Descartes_variations']==2 and roots['positive_root_count']==2
        assert roots['completeness_proof']=='Descartes_upper_equals_exact_disjoint_bracket_lower'
        assert roots['Sturm_count'] is None
        sq=QP(S);oldz=Q(0);oldt=Q(0);verified=[]
        locked_a=[Q(62145027,2**32),Q(92903691,2**31)]
        locked_b=[Q(15536257,2**30),Q(185807383,2**32)]
        for i,row in enumerate(roots['roots']):
            a,b=Q(row['z_lower']),Q(row['z_upper'])
            ta,tb=Q(row['T_lower']),Q(row['T_upper'])
            assert 0<oldz<a<b and oldt<ta<tb
            assert sq(a)*sq(b)<0
            owners=[(k,m) for k,(f,m) in enumerate(sf) if QP(f)(a)*QP(f)(b)<0]
            assert owners==[(row['squarefree_factor_index'],row['multiplicity'])]
            assert row['multiplicity']==1 and ta==locked_a[i] and tb==locked_b[i]
            assert ta*ta<=a<b<=tb*tb
            oldz=b;oldt=tb
            verified.append({'T_lower':str(ta),'T_upper':str(tb),'z_lower':str(a),'z_upper':str(b),
                             'multiplicity':1,'opposite_exact_signs':True})
        assert len(verified)==2 and roots['positive_roots_with_multiplicity']==2
        # Fresh independent pi enclosure; evaluating a frozen Gamma is not a full G call.
        ctx.prec=768;pl,pu=endpoints(arb.pi())
        assert F(str(ell))<pl<pu<F(str(u))
        def saved_pi_endpoint(which):
            row=constants['pi'];return F(int(row[which+'_mantissa']))*F(2)**int(row[which+'_binary_exponent'])
        spl,spu=saved_pi_endpoint('lower'),saved_pi_endpoint('upper')
        assert F(str(ell))<spl<=pl<=pu<=spu<F(str(u))
        gamma=F(f23['Gamma_G']);assert gamma>0
        Glo,Ghi=gamma/(128*pu**5),gamma/(128*pl**5)
        result={'status':'PASS_INDEPENDENT_LOCAL_N_D_P_IDENTITY_AND_TWO_COMPLETE_POSITIVE_ROOTS',
                'local_principal_slots':417,'positive_poles':359,'N_degree':414,'D_degree':417,
                'D_positive':True,'N_D_coprime':True,'P_positive_scale_verified':True,
                'dense_D_rebuilt':False,'new_global_numerator':False,'Descartes_variations':2,
                'positive_roots':verified,'squarefree_verified':True,'root_isolation_calls':0,'Sturm_calls':0,
                'ell':str(ell),'upper_pi':str(u),'c':str(c),'pi_768_exact_bounds':[str(pl),str(pu)],
                'G_from_F023_exact_bounds':[str(Glo),str(Ghi)],'new_G_full_calls':0,
                'input_sha256':INPUTS,'script_sha256':sha(Path(__file__)),
                'pid':os.getpid(),'elapsed_seconds':time.time()-start}
        write('root_identity_results.json',result)
        receipt['status']='ROOT_IDENTITY_PASS_NATIVE_NOT_STARTED'
        receipt['root_identity_seconds']=time.time()-start
    except BaseException as e:
        receipt['status']='FAILED_STOPPED';receipt['exception']=repr(e);raise
    finally:write('audit_receipt.json',receipt);signal.alarm(0)
    print(json.dumps({k:result[k] for k in ('status','pid','elapsed_seconds','local_principal_slots','dense_D_rebuilt')}))
if __name__=='__main__':main()
