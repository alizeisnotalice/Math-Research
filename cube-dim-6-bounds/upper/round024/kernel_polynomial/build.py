#!/usr/bin/env python3
"""One explicitly released exact N/D/P construction/isolation batch."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import argparse,hashlib,inspect,json,math,os,resource,subprocess,time
from fractions import Fraction as F
from pathlib import Path
from flint import fmpq as Q,fmpq_poly as QP,fmpz_poly as ZP,arb,ctx
import flint,sympy
from sympy import Poly,Symbol,ZZ,Rational,oo
from sympy.polys import rootisolation,densearith,densetools,polytools

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[2]/'upper/round023/partial_fraction_kernel'
AUDIT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/round023/partial_fraction_review')
F_REPORT='2851b03ea55e4ce3e95e7bbbae0b97daa9724a01f0de051aadc586255427544e'
F_MANIFEST='c8b3f27cb4a4cada0993e50e89366b7a5f9a4997946f3986a70a31093cf9d70a'
B_TABLE='8446b3fe84136c005d66501a64a7978e1c1adc452bc1993acd25545ddc8ff299'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path,obj):
    path=Path(path);temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n');os.replace(temp,path)

def encode(poly):
    return {'integer_coefficients_ascending':[str(v)for v in poly.numer().coeffs()],
            'positive_common_denominator':str(poly.denom())}

def decode(obj):
    d=int(obj['positive_common_denominator']);assert d>0
    return QP([Q(int(v),d)for v in obj['integer_coefficients_ascending']])

def zdecode(coeffs):return ZP([int(v)for v in coeffs])

def interval(value):
    a,e=value.lower().man_exp();b,f=value.upper().man_exp()
    return {'lower_mantissa':str(a),'lower_binary_exponent':str(e),
            'upper_mantissa':str(b),'upper_binary_exponent':str(f)}

def bounds(obj):
    return (F(int(obj['lower_mantissa']))*F(2)**int(obj['lower_binary_exponent']),
            F(int(obj['upper_mantissa']))*F(2)**int(obj['upper_binary_exponent']))

def variations(coeffs):
    signs=[1 if a>0 else -1 for a in coeffs if a]
    return sum(a!=b for a,b in zip(signs,signs[1:]))

def to_sympy(poly):
    return Poly.from_list([int(a)for a in reversed(poly.coeffs())],Symbol('z'),domain=ZZ)

def qeval(poly,x):return QP([Q(int(v))for v in poly.coeffs()])(Q(x.numerator,x.denominator))

def register():
    assert not(HERE/'input_manifest.json').exists()
    assert sha(BASE/'partial_tables.json')==B_TABLE
    assert sha(AUDIT/'report.md')==F_REPORT and sha(AUDIT/'manifest.json')==F_MANIFEST
    files=[BASE/f for f in ('report.md','manifest.json','partial_tables.json','results.json','theorem.md','checked_results.json')]
    files += [AUDIT/'report.md',AUDIT/'manifest.json',HERE/'registration.md',HERE/'build.py',HERE/'check.py']
    files += [Path(inspect.getfile(m))for m in (rootisolation,densearith,densetools,polytools)]
    flint_home=Path(flint.__file__).parent
    files += list((flint_home/'.dylibs').glob('libflint*.dylib'))
    files += list((flint_home/'.dylibs').glob('libgmp*.dylib'))
    files += list((flint_home/'.dylibs').glob('libmpfr*.dylib'))
    obj={'status':'REGISTERED_AFTER_EXPLICIT_MAIN_RELEASE_AND_F023_PASS','pid':os.getpid(),
         'batch_count_cap':1,'seconds_cap':1200,'degree_cap':512,'positive_bracket_cap':4,
         'sampled_rss_cap':768*1024**2,'new_K_H_native_calls':0,'new_G_full_calls':0,
         'old_G_Arb_enclosure_cap':1,'pi_Arb_precision':512,'root_eps':'1/1208925819614629174706176',
         'T_bracket_bit_schedule':[32,64,128,256,512],
         'fixed_c':'1/4100','pi_lower':'3141592653/1000000000','pi_upper':'1570796327/500000000',
         'F023_report_sha256':F_REPORT,'F023_manifest_sha256':F_MANIFEST,
         'inputs':[{'path':str(p),'sha256':sha(p)}for p in files],
         'dependencies':{'python':sys.version,'python_flint':flint.__version__,'sympy':sympy.__version__},
         'API_evidence':{'count_roots':'SymPy dup_count_real_roots builds exact Sturm chain and counts distinct roots',
                         'intervals':'SymPy exact VAS rational isolation; no approximate roots'}}
    dump(HERE/'input_manifest.json',obj)
    print(json.dumps({'status':obj['status'],'input_manifest_sha256':sha(HERE/'input_manifest.json'),'pid':os.getpid()},indent=2))

def worker(start_epoch):
    manifest=json.loads((HERE/'input_manifest.json').read_text())
    for r in manifest['inputs']:assert sha(r['path'])==r['sha256']
    history=[]
    def checkpoint(stage,status,extra=None):
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
        row={'stage':stage,'status':status,'elapsed_seconds':time.time()-start_epoch,'rss_bytes':rss}
        if extra:row.update(extra)
        history.append(row);dump(HERE/'checkpoint.json',{'worker_pid':os.getpid(),'stages':history,'last_stage':row})
        print(json.dumps(row),flush=True)
        assert row['elapsed_seconds']<manifest['seconds_cap'],'batch time gate'
        assert rss<manifest['sampled_rss_cap'],'resident memory gate'
    try:
        checkpoint('constants','started')
        ctx.prec=512
        ell=Q(3141592653,10**9);up=Q(3141592654,10**9);c=Q(1,4100)
        pi=arb.pi();pi_bounds=interval(pi);pl,pu=bounds(pi_bounds)
        assert F(str(ell))<pl<=pu<F(str(up))
        old=json.loads((BASE/'results.json').read_text());gamma=Q(old['Gamma_G'])
        G=arb(gamma)/(128*pi**5);gb=interval(G);gl,gu=bounds(gb)
        assert gl>0
        safe_g=(F(str(gamma))/(128*pu**5),F(str(gamma))/(128*pl**5))
        tail=F(math.ceil(max(gu,safe_g[1])*4100*10**12),10**12)
        constants={'precision_bits':512,'pi':pi_bounds,'ell':str(ell),'upper_pi':str(up),'c':str(c),
                   'strict_pi_outer_bounds_verified':True,'G_from_old_Gamma':gb,'Gamma_source_sha256':sha(BASE/'results.json'),
                   'safe_G_rational_bounds_from_pi':[str(v)for v in safe_g],
                   'large_T_tail_bound':str(tail),'tail_inequality':'T>=tail => K_T/T<=G/T<=c',
                   'old_G_enclosures':1,'new_G_full_calls':0,'new_kernel_calls':0}
        dump(HERE/'constants.json',constants);checkpoint('constants','complete')
        checkpoint('synthesis','started')
        table=json.loads((BASE/'partial_tables.json').read_text())
        factors={Q(r['b_squared']):r['exponent']for r in table['denominator_factors']}
        assert sum(factors.values())==417 and all(q>0 for q in factors)
        local={q:{}for q in factors}
        for r in table['real_z_principal']:
            local[Q(r['b_squared'])][r['order']]=Q(r['coefficient'])
        pairs=[]
        for q,m in sorted(factors.items()):
            w=QP([q,1]);n=QP([0])
            for j,a in local[q].items():n+=a*w**(m-j)
            pairs.append((n,w**m))
        level=0
        while len(pairs)>1:
            out=[]
            for i in range(0,len(pairs),2):
                if i+1==len(pairs):out.append(pairs[i]);continue
                n1,d1=pairs[i];n2,d2=pairs[i+1]
                out.append((n1*d2+n2*d1,d1*d2))
            pairs=out;level+=1
            checkpoint('synthesis','balanced_level',{'level':level,'pairs':len(pairs)})
        N,D=pairs[0]
        assert D.degree()==417<=manifest['degree_cap'] and D.leading_coefficient()==1
        assert N[0]==0 and N.degree()<=414 and N.gcd(D).degree()==0
        assert N.leading_coefficient()==Q(old['full_weight_checksum']) and N.degree()==414
        kernel={'N':encode(N),'D_factors':[{'s':str(q),'exponent':m}for q,m in sorted(factors.items())],
                'D_degree':D.degree(),'N_degree':N.degree(),'N_zero_at_zero':True,'gcd_N_D_one':True,
                'D_positive_on_positive_axis':True,'partial_table_sha256':B_TABLE,
                'construction':'balanced exact local rational sum','D_coefficient_binding_sha256':hashlib.sha256(json.dumps(encode(D),sort_keys=True).encode()).hexdigest()}
        dump(HERE/'kernel_polynomials.json',kernel);checkpoint('synthesis','complete',{'degree_N':N.degree(),'degree_D':D.degree()})
        checkpoint('critical_squarefree','started')
        P=N-64*ell**6*c*D
        den=int(P.denom());numer=P.numer();content=int(numer.content());assert den>0 and content>0
        zpoly=ZP([int(a)//content for a in numer.coeffs()])
        scale=Q(content,den);assert QP(zpoly)*scale==P
        assert zpoly.degree()==417<=manifest['degree_cap'] and zpoly[0]<0 and zpoly.leading_coefficient()<0
        sf_content,sf=zpoly.factor_squarefree();rebuild=ZP([sf_content]);S=ZP([1])
        for f,m in sf:
            assert f.gcd(f.derivative()).degree()==0
            rebuild*=f**m;S*=f
        assert rebuild==zpoly and S[0]!=0 and S.gcd(S.derivative()).degree()==0
        critical={'primitive_integer_P_ascending':[str(a)for a in zpoly.coeffs()],
                  'positive_scale_P_equals_scale_times_primitive':str(scale),'ell':str(ell),'c':str(c),
                  'squarefree_content':str(sf_content),'squarefree_factors':[{'multiplicity':m,'integer_coefficients_ascending':[str(a)for a in f.coeffs()]}for f,m in sf],
                  'squarefree_part_ascending':[str(a)for a in S.coeffs()],
                  'P_zero_value_negative':True,'P_leading_negative':True,
                  'kernel_polynomials_sha256':sha(HERE/'kernel_polynomials.json')}
        dump(HERE/'critical_polynomial.json',critical);checkpoint('critical_squarefree','complete',{'P_degree':zpoly.degree(),'S_degree':S.degree(),'sqf_factor_count':len(sf)})
        checkpoint('positive_roots','started')
        V=variations(list(reversed(S.coeffs())));SS=to_sympy(S);sturm_count=None
        if V>4:
            sturm_count=int(SS.count_roots(0,oo))
            checkpoint('positive_roots','exact_sturm_count',{'count':sturm_count,'Descartes_upper':V})
            assert sturm_count<=4,'positive root count gate'
        found=SS.intervals(inf=0,eps=Rational(1,2**80),sqf=True)
        assert len(found)<=4,'positive bracket gate'
        lower_bound=F(1,1)/(1+max(F(abs(int(a)),abs(int(S[0])))for a in S.coeffs()[1:]))
        brackets=[]
        for a,b in found:
            a=F(int(a.p),int(a.q));b=F(int(b.p),int(b.q))
            assert 0<=a<=b and b>0
            if a==0:a=lower_bound/2
            assert a>0 and a<=b
            if a==b:assert qeval(S,a)==0
            else:assert qeval(S,a)*qeval(S,b)<0
            brackets.append((a,b))
        assert all(b<=aa for (a,b),(aa,bb)in zip(brackets,brackets[1:]))
        if sturm_count is None and len(brackets)==V:
            proof='Descartes_upper_equals_exact_disjoint_bracket_lower'
            count=len(brackets)
        else:
            if sturm_count is None:sturm_count=int(SS.count_roots(0,oo))
            count=sturm_count;proof='exact_Sturm_count_roots'
            assert count<=4 and count==len(brackets)
        roots=[]
        for a,b in brackets:
            owners=[]
            for index,(f,m)in enumerate(sf):
                if (qeval(f,a)==0 if a==b else qeval(f,a)*qeval(f,b)<0):owners.append((index,m))
            assert len(owners)==1
            roots.append({'z_lower':str(a),'z_upper':str(b),'squarefree_factor_index':owners[0][0],'multiplicity':owners[0][1]})
        for bits in (32,64,128,256,512):
            Sgrid=1<<bits;converted=[]
            for root in roots:
                a,b=F(root['z_lower']),F(root['z_upper'])
                an,ad=math.isqrt(a.numerator),math.isqrt(a.denominator)
                if a==b and an*an==a.numerator and ad*ad==a.denominator:
                    lo=hi=F(an,ad)
                else:
                    lo=F(math.isqrt(a.numerator*Sgrid*Sgrid//a.denominator),Sgrid)
                    if lo==0:lo=a/(a+1)
                    hi=F(math.isqrt(b.numerator*Sgrid*Sgrid//b.denominator)+1,Sgrid)
                assert 0<lo<=hi and lo*lo<=a<=b<=hi*hi
                converted.append((lo,hi))
            if all(v<=w for (u,v),(w,x)in zip(converted,converted[1:])):break
        else:raise AssertionError('T bracket bit gate')
        for root,(lo,hi)in zip(roots,converted):root.update({'T_lower':str(lo),'T_upper':str(hi),'T_bracket_bits':bits})
        root_data={'status':'EXACT_ALL_POSITIVE_P_ROOTS_ISOLATED','positive_root_count':count,
                   'positive_roots_with_multiplicity':sum(r['multiplicity']for r in roots),'roots':roots,
                   'Descartes_variations':V,'Sturm_count':sturm_count,'completeness_proof':proof,
                   'S_zero_nonzero':True,'reciprocal_Cauchy_nonzero_root_lower_bound':str(lower_bound),
                   'critical_polynomial_sha256':sha(HERE/'critical_polynomial.json'),
                   'no_positive_root_consequence':('P<0 on z>0, hence p_T<c and K_T/T<c; fixed-model obstruction pending full source audit'if count==0 else None),
                   'K_evaluations_at_roots':0,'all_time_claim_when_roots_exist':False}
        dump(HERE/'roots.json',root_data);checkpoint('positive_roots','complete',{'count':count,'proof':proof})
        checkpoint('saved_checker','started')
        import check
        checked=check.saved(start_epoch+manifest['seconds_cap'])
        checkpoint('saved_checker','complete',{'status':checked['status']})
        dump(HERE/'results.json',{'status':'B024_EXACT_CONSTRUCTION_ISOLATION_CHECKED_STOPPED','degree_N':N.degree(),'degree_D':D.degree(),
             'positive_root_count':count,'positive_roots_with_multiplicity':root_data['positive_roots_with_multiplicity'],
             'root_completeness_proof':proof,'new_kernel_calls':0,'new_G_full_calls':0,'old_G_enclosures':1,
             'new_C6_bound':False,'root_K_values_computed':False,'elapsed_seconds':time.time()-start_epoch})
    except BaseException as exc:
        dump(HERE/'failure.json',{'status':'B024_FAILED_STOPPED','error':repr(exc),'last_stage':history[-1]if history else None,
                                'no_automatic_retry':True,'new_kernel_calls':0})
        raise

def run():
    assert not(HERE/'run_receipt.json').exists(),'one registered batch; no automatic replay'
    manifest=json.loads((HERE/'input_manifest.json').read_text())
    for row in manifest['inputs']:assert sha(row['path'])==row['sha256']
    start=time.time();record={'status':'BATCH_STARTED_AFTER_EXPLICIT_RELEASE','driver_pid':os.getpid(),'started_epoch':start,
        'batch_count':1,'new_kernel_calls':0,'input_manifest_sha256':sha(HERE/'input_manifest.json')}
    dump(HERE/'run_receipt.json',record)
    command=[sys.executable,str(HERE/'build.py'),'worker','--start',str(start)]
    process=subprocess.Popen(command,cwd=HERE)
    record['worker_pid']=process.pid;dump(HERE/'run_receipt.json',record)
    try:
        code=process.wait(timeout=manifest['seconds_cap'])
    except subprocess.TimeoutExpired:
        process.kill();code=process.wait();record['timeout']=True
    except BaseException:
        if process.poll()is None:process.kill()
        code=process.wait();record['interrupted']=True
        raise
    finally:
        if process.poll()is None:process.kill();process.wait()
        record['returncode']=process.returncode;record['elapsed_seconds']=time.time()-start
        record['status']='BATCH_EXIT0_STOPPED'if process.returncode==0 else 'BATCH_FAILED_OR_TIMEOUT_STOPPED'
        record['all_children_waited']=True;dump(HERE/'run_receipt.json',record)
    if code:raise SystemExit(2)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('register','run','worker'));p.add_argument('--start',type=float)
    args=p.parse_args()
    if args.mode=='register':register()
    elif args.mode=='run':run()
    else:worker(args.start)
