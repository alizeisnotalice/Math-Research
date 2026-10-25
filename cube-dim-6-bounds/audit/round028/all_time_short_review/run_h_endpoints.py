#!/usr/bin/env python3
"""Exactly two locked C010 H512 calls, independently paying critical brackets."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,subprocess,sys,time
sys.set_int_max_str_digits(0)

HERE=Path(__file__).resolve().parent
TOP=Path('/Users/zhengzhihao/Desktop/T')
CROOT=TOP/'Math-Research-cube6-audit/cube-dim-6-bounds/upper'
ENGINE=CROOT/'round010/arb_engine'
C26=CROOT/'round026/integrated_pf_potential'
FAN=TOP/'Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round003/direct_kernel/fan_integer_input.txt'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def write(n,d):(HERE/n).write_text(json.dumps(d,indent=2)+'\n')
def walk(obj,path):
    for key in path:obj=obj[key]
    return obj
def bounds(d):
    if 'common_binary_exponent' in d:
        factor=F(2)**int(d['common_binary_exponent'])
        return int(d['lower_mantissa'])*factor,int(d['upper_mantissa'])*factor
    return tuple(F(int(d[s+'_mantissa']))*F(2)**int(d[s+'_binary_exponent']) for s in ('lower','upper'))

def main():
    locked=load(HERE/'locked_inputs.json');roots=load(HERE/'root_identity_results.json')
    receipt=load(HERE/'audit_receipt.json')
    assert receipt['status']=='ROOT_IDENTITY_PASS_NATIVE_NOT_STARTED'
    assert receipt['native_H_attempts']==0
    assert not (HERE/'native_endpoint_results.json').exists()
    assert time.time()<receipt['deadline_unix']
    for name,h in locked['review_sha256'].items():assert sha(HERE/name)==h,name
    for path,h in locked['frozen_sha256'].items():assert sha(path)==h,path
    for path,h in roots['input_sha256'].items():assert sha(path)==h,path
    c26=load(C26/'manifest.json')
    for path,h in c26['external_readonly_sha256'].items():assert sha(path)==h,path
    assert sha(C26/'pf_integral.py')=='5f0fe640bfe73fca40fcad55fcc55bdb6b066a5724664f1fa23abb62302013ad'
    deps=load(ENGINE/'dependencies.json')
    for path,h in deps['runtime_library_sha256'].items():assert sha(path)==h,path
    binary=ENGINE/'scratch/arb_kernel'
    assert sha(binary)=='429cdccca45a69c1de6c04725231297e00ee8d3951ccbb1a6141feb934f13bb8'
    assert sha(ENGINE/'arb_kernel.cpp')=='e857e480efae68ec99fd82f97c64f4756b6bff465476e3a788d9c283f5b751b7'
    assert sha(FAN)=='7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'
    # Exact safe bounds for every permuted/signed numerator; no kernel traversal.
    words=iter(map(int,FAN.read_text().split()));assert next(words)==21
    max_n=max_d=0
    for _ in range(21):
        assert next(words)>0
        for _ in range(6):
            q=next(words);v=[next(words) for _ in range(6)];assert q>0
            max_n=max(max_n,50*sum(map(abs,v)));max_d=max(max_d,50*q)
    assert next(words,None) is None
    assert max_n<2**63 and max_d<2**63 and max_n*max_d<2**127
    expected_b=[F(15536257,2**30),F(185807383,2**32)]
    assert len(roots['positive_roots'])==2 and len(locked['PF_crosschecks'])==2
    Glo,Ghi=map(F,roots['G_from_F023_exact_bounds'])
    ell,u,c=F(roots['ell']),F(roots['upper_pi']),F(roots['c'])
    assert c==F(1,4100);alpha_upper=(u/ell)**6
    inputs={'binary_sha256':sha(binary),'source_sha256':sha(ENGINE/'arb_kernel.cpp'),
            'fan_sha256':sha(FAN),'runtime_sha256':deps['runtime_library_sha256'],
            'root_identity_sha256':sha(HERE/'root_identity_results.json'),
            'independent_G_source':'F023 exact Gamma and fresh 768-bit pi in root audit',
            'max_raw_node_numerator_bound':max_n,'max_raw_node_denominator':max_d,
            'script_sha256':sha(Path(__file__)),
            'calls':[{'T_upper':str(t),'timeout_seconds':120} for t in expected_b]}
    write('native_input_manifest.json',inputs)
    records=[]
    try:
        for i,(root,pcfg) in enumerate(zip(roots['positive_roots'],locked['PF_crosschecks'])):
            assert time.time()<receipt['deadline_unix']
            a,b=F(root['T_lower']),F(root['T_upper']);assert b==expected_b[i]
            pfraw=load(pcfg['file']);pf=walk(pfraw,pcfg['record_path'])
            assert F(walk(pf,pcfg['T_path']))==b
            pl,pu=bounds(walk(pf,pcfg['K_path']));assert 0<pl<=pu
            cmd=[str(binary),'full',str(FAN),str(b.numerator),str(b.denominator),
                 '50','50','50','50','45','27','0']
            receipt['status']='NATIVE_ENDPOINT_AUDIT_RUNNING'
            receipt['native_H_attempts']+=1
            receipt.setdefault('commands',[]).append(cmd)
            write('audit_receipt.json',receipt)
            tick=time.time();env=dict(os.environ)
            for key in ('DYLD_LIBRARY_PATH','DYLD_FALLBACK_LIBRARY_PATH','DYLD_INSERT_LIBRARIES'):
                env.pop(key,None)
            run=subprocess.run(cmd,capture_output=True,text=True,timeout=120,env=env)
            wall=time.time()-tick
            (HERE/f'H_root{i+1}.json').write_text(run.stdout)
            receipt.setdefault('native_exit_codes',[]).append(run.returncode)
            receipt.setdefault('native_wall_seconds',[]).append(wall)
            write('audit_receipt.json',receipt)
            if run.returncode or run.stderr:
                (HERE/f'H_root{i+1}.stderr.txt').write_text(run.stderr)
                raise RuntimeError('Native returned error or stderr; stop without retry')
            native=json.loads(run.stdout)
            assert native['flint_version']=='3.6.0' and native['precision_bits']==512
            assert native['full_group_terms']==967680 and native['computed_DD_evaluations']==483840
            assert native['zero_node_terms']==88320 and native['confluent_terms']==201408
            assert native['exact_global_sign_pairing'] is True
            hlo,hhi=bounds(native['dyadic_interval']);assert 0<hlo<=hhi
            klo,khi=Glo-hhi,Ghi-hlo;assert 0<klo<=khi
            margin=c*a-alpha_upper*khi
            assert margin>0,'Independent root payment failed; stop without retry'
            assert max(klo,pl)<=min(khi,pu),'PF/reference crosscheck disjoint; stop'
            records.append({'root_index':i+1,'a':str(a),'b':str(b),
                            'H_exact_bounds':[str(hlo),str(hhi)],
                            'K_independent_exact_bounds':[str(klo),str(khi)],
                            'PF_exact_bounds':[str(pl),str(pu)],
                            'PF_strictly_inside_independent_K':klo<pl<=pu<khi,
                            'critical_F_lower_exact':str(margin),'critical_F_lower_display_only':float(margin),
                            'strict_payment':True,'native_wall_seconds':wall,
                            'native_sha256':sha(HERE/f'H_root{i+1}.json')})
            write('native_endpoint_results.json',{'status':'PARTIAL_NOT_YET_ALL_TIME',
                                                 'native_H_calls':receipt['native_H_attempts'],'roots':records})
        assert len(records)==receipt['native_H_attempts']==2
        out={'status':'PASS_TWO_INDEPENDENT_G_MINUS_H_ROOT_PAYMENTS',
             'native_H_calls':2,'new_PF_full_sums':0,'new_G_full_calls':0,'Sturm_calls':0,
             'roots':records,'alpha_upper_exact':str(alpha_upper),'c':'1/4100',
             'supported_all_positive_times':True,'unrestricted_probability_time_mixtures':True,
             'model_cost_lower':'1025/16','new_C6_bound':False,
             'receipt_start_unix':receipt['start_unix'],'elapsed_since_audit_start':time.time()-receipt['start_unix']}
        write('native_endpoint_results.json',out)
        receipt['status']='INDEPENDENT_ALL_TIME_AUDIT_PASS_STOPPED'
    except BaseException as e:
        receipt['status']='FAILED_STOPPED_NO_RETRY';receipt['exception']=repr(e);raise
    finally:
        receipt['elapsed_seconds']=time.time()-receipt['start_unix'];write('audit_receipt.json',receipt)
    print(json.dumps({'status':out['status'],'native_H_calls':2,
                      'margins_display_only':[r['critical_F_lower_display_only'] for r in records],
                      'native_wall_seconds':receipt['native_wall_seconds'],'total_audit_seconds':receipt['elapsed_seconds']}))
if __name__=='__main__':main()
