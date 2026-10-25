#!/usr/bin/env python3
"""Registered beta=1, T=2/45 at at most two existing exact points."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as Q
from time import perf_counter
import ast, hashlib, itertools, json, math, os, signal, subprocess
import flint
from flint import arb, ctx, fmpq, fmpz
from audit_frame import HERE, PROJECT, NEW_FAN, B010, ENGINE, BINARY, H_SOURCE, T, sha, load_module

ctx.prec=512
def A(q):return arb(q.numerator)/q.denominator
def ball(v):return {'ball':v.str(130,more=True),'lower':v.lower().str(130,more=True),
                   'upper':v.upper().str(130,more=True)}
def native_ball(row):
    p=row['dyadic_interval'];e=int(p['common_binary_exponent'])
    lo=int(p['lower_mantissa']);hi=int(p['upper_mantissa']);assert lo<=hi
    return arb((lo+hi,e-1),(hi-lo,e-1))
def write(name,value):(HERE/name).write_text(json.dumps(value,indent=2)+'\n')

def main():
    start=perf_counter()
    manifest=json.loads((HERE/'input_manifest.json').read_text())
    for row in manifest['inputs']:assert sha(row['path'])==row['sha256'],row['path']
    audit=json.loads((HERE/'fan_audit.json').read_text())
    assert audit['status']=='PASS_EXACT_AUGMENTED_FRAME_AND_FAN'
    assert audit['new_fan_sha256']==sha(NEW_FAN)
    assert audit['manifest_sha256']==sha(HERE/'input_manifest.json')
    inputs=[HERE/'registration.md',HERE/'theorem.md',HERE/'input_manifest.json',HERE/'fan_audit.json',
            HERE/'verify_exact_fan_results.json',NEW_FAN,Path(__file__),HERE/'audit_frame.py']
    kernel_manifest={'status':'REGISTERED_BEFORE_FIRST_NEW_FRAME_PASS','beta':'1','T':str(T),
        'max_full_passes':4,'max_seconds':1200,'points':manifest['points'],
        'inputs':[{'path':str(p),'sha256':sha(p)}for p in inputs]}
    write('kernel_manifest.json',kernel_manifest)
    out={'status':'RUNNING_REGISTERED_AUGMENTED_CANDIDATE','pid':os.getpid(),'beta':'1','T':str(T),
         'points':[],'new_G_passes':0,'new_H_passes':0,'new_full_passes':0,
         'python':sys.version,'python_flint':flint.__version__,'arb_bits':512,
         'kernel_manifest_sha256':sha(HERE/'kernel_manifest.json'),
         'source_sha256':sha(__file__),'new_C6_bound':False,'global_m_certified':False}
    def save():
        out['elapsed_seconds']=perf_counter()-start;write('results.json',out)
    signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(TimeoutError('registered 1200-second kernel-stage cap')))
    signal.alarm(1200)
    try:
        # Import only audited function definitions, never the B010 executable top level.
        source=B010/'evaluate_point.py';tree=ast.parse(source.read_text())
        names={'generic_dd','laurent_dd','anchors','parse_fan','exact_green'}
        definitions=[node for node in tree.body if isinstance(node,ast.FunctionDef)and node.name in names]
        assert {node.name for node in definitions}==names
        env={'fmpq':fmpq,'arb':arb,'math':math,'itertools':itertools,'perf_counter':perf_counter,
             'ball':ball,'FAN':NEW_FAN,'__builtins__':__builtins__}
        exec(compile(ast.Module(body=definitions,type_ignores=[]),str(source),'exec'),env)
        out['G_function_source_sha256']=sha(source)
        out['G_functions_loaded_without_top_level_execution']=sorted(names)
        out['scalar_confluent_anchors']=env['anchors']();save()
        gcones=env['parse_fan']()
        for point in manifest['points']:
            index=point['point_index'];x=point['numerators'];den=point['denominator']
            assert x[0]==den and all(x[i]>x[i+1]>0 for i in range(5))
            row={'point_index':index,'point_numerators':x,'denominator':den,
                 'off_new_S2_by_six_distinct_positive_coordinates':True,'new_full_passes':0}
            out['points'].append(row);save()
            g=env['exact_green'](x,den,gcones)
            g['G_arb512']=g.pop('G_arb384');g['actual_precision_bits']=512
            row['G']=g;row['new_full_passes']+=1;out['new_G_passes']+=1;out['new_full_passes']+=1;save()
            G=arb(fmpq(fmpz(g['G_times_pi5_numerator']),fmpz(g['G_times_pi5_denominator'])))/arb.pi()**5
            assert (G-arb(g['G_arb512']['ball'])).contains(0)
            print(json.dumps({'point_index':index,'new_G':g['G_arb512']['ball'],'seconds':g['elapsed_seconds']}),flush=True)
            if index==1:
                safety=audit['native_input_range_checks'][0]
                assert safety['native_int64_safe']and safety['native_cross_product_int128_safe']
                args=[str(BINARY),'full',str(NEW_FAN),str(T.numerator),str(T.denominator),str(den),*map(str,x)]
                tic=perf_counter();proc=subprocess.run(args,text=True,capture_output=True,
                    timeout=max(1,1200-(perf_counter()-start)),check=True)
                raw=json.loads(proc.stdout);write('point1_native_H.json',raw)
                assert raw['precision_bits']==512 and raw['full_group_terms']==967680
                H=native_ball(raw)
                row['H']={'method':'C010 native Arb512 with exact dyadic output','native_file':'point1_native_H.json',
                          'native_sha256':sha(HERE/'point1_native_H.json'),'H_arb512':ball(H),
                          'command':args,'elapsed_seconds':perf_counter()-tic,
                          'complete_group_terms':raw['full_group_terms'],'zero_node_terms':raw['zero_node_terms'],
                          'confluent_terms':raw['confluent_terms'],'native_returncode':proc.returncode}
            else:
                assert not audit['native_input_range_checks'][1]['native_int64_safe']
                m=load_module(H_SOURCE,'c012_arbitrary_integer_H');m.FAN=NEW_FAN;ctx.prec=512
                hr=m.evaluate(x,den,A(T),m.fan(),start,gradient=False)
                H=arb(hr['H_arb384'])
                row['H']={'method':'existing Python arbitrary-integer rational nodes and Arb at 512 bits',
                          'actual_precision_bits':512,'source_sha256':sha(H_SOURCE),'H_arb512':ball(H),
                          'original_serialized_H':hr['H_arb384'],'native_engine_called':False,
                          'complete_group_terms':hr['complete_group_terms'],
                          'computed_sign_paired_terms':hr['computed_sign_paired_terms'],
                          'elapsed_seconds':hr['elapsed_seconds']}
            assert row['H']['complete_group_terms']==967680
            row['new_full_passes']+=1;out['new_H_passes']+=1;out['new_full_passes']+=1
            K=G-H;assert H>0 and K>0
            threshold=A(T/4096);m0=A(T/4000)
            row.update({'K_beta_T':ball(K),'necessary_threshold':str(T/4096),'sufficient_62_5_target':str(T/4000),
                        'K_minus_T_over_4096':ball(K-threshold),'K_minus_T_over_4000':ball(K-m0),
                        'necessary_condition_for_below_64_passes':bool(K>threshold),
                        'rules_out_below_64_candidate':bool(K<=threshold),
                        'rules_out_62_5_sufficient_target':bool(K<m0),
                        'formal_point_cost':ball(A(T)/(64*K))})
            save();print(json.dumps({'point_index':index,'K':row['K_beta_T'],'cost':row['formal_point_cost'],
                'rules_out_below_64_candidate':row['rules_out_below_64_candidate'],'new_full_passes':out['new_full_passes']},indent=2),flush=True)
            if K<=threshold:
                out['status']='STOP_FIRST_CERTIFIED_NECESSARY_CONDITION_FAILURE';out['stopped_at_point']=index;break
            if not K>threshold:
                out['status']='STOP_UNRESOLVED_THRESHOLD_ENCLOSURE';out['stopped_at_point']=index;break
        else:out['status']='COMPLETE_TWO_POINTS_NECESSARY_CONDITIONS_ONLY'
        assert out['new_full_passes']<=4
        out['exit_code']=0;out['no_more_points_beta_or_times_permitted']=True;save()
    except BaseException as ex:
        out.update({'status':'ERROR_STOPPED','error':repr(ex),'exit_code':1});save();raise
    finally:signal.alarm(0)

if __name__=='__main__':main()
