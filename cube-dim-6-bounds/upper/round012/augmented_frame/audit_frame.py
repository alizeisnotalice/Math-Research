#!/usr/bin/env python3
"""Exact new rates on the old rational fan; no kernel evaluations."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from time import perf_counter
import hashlib, importlib.util, json, math, os, signal
from flint import fmpq_mat

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
OLD_FAN=PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
GEOMETRY=PROJECT/'upper/round002/fan/exact_cells_results.json'
VERIFY=PROJECT/'upper/round002/fan/verify_exact_fan.py'
OLD_VERIFY=PROJECT/'upper/round002/fan/verify_exact_fan_results.json'
OLD_GENERATE=PROJECT/'upper/round003/direct_kernel/generate_integer_input.py'
B010=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round010/short_minimum_diagnostic')
ENGINE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round010/arb_engine')
BINARY=ENGINE/'scratch/arb_kernel'
H_SOURCE=PROJECT/'upper/round007/g_time_moments/run_moments.py'
NEW_FAN=HERE/'fan_beta1.txt'
T=Q(2,45)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,value):(HERE/name).write_text(json.dumps(value,indent=2)+'\n')
def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def parse_fan(path):
    it=iter(Path(path).read_text().split());cones=[]
    for _ in range(int(next(it))):
        determinant=int(next(it));rays=[]
        for _ in range(6):
            q=int(next(it));rays.append((q,tuple(int(next(it))for _ in range(6))))
        cones.append((determinant,rays))
    assert next(it,None) is None
    return cones

def main():
    start=perf_counter()
    paths=[HERE/'registration.md',HERE/'theorem.md',Path(__file__),HERE/'evaluate.py',
        OLD_FAN,GEOMETRY,VERIFY,OLD_VERIFY,OLD_GENERATE,H_SOURCE,
        PROJECT/'upper/round002/kernel/002_kernel_interface.md',
        B010/'evaluate_point.py',B010/'point1.json',B010/'point2.json',
        B010/'point1_results.json',B010/'other_time_results.json',B010/'batch2_results.json',
        ENGINE/'arb_kernel.cpp',ENGINE/'build_record.json',ENGINE/'validation_results_v3.json',
        ENGINE/'dependencies.json',BINARY]
    points=[]
    for i in (1,2):
        row=json.loads((B010/f'point{i}.json').read_text())
        points.append({'point_index':i,'denominator':row['denominator'],
                       'numerators':row['candidate_numerators']})
    manifest={'status':'REGISTERED_BEFORE_FRAME_AUDIT_AND_NEW_KERNELS','beta':'1','T':str(T),
        'points':points,'max_new_full_passes':4,'max_seconds_per_stage':1200,
        'inputs':[{'path':str(p),'sha256':sha(p)}for p in paths]}
    write('input_manifest.json',manifest)
    signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(TimeoutError('1200-second frame audit')))
    signal.alarm(1200)
    try:
        verify=load_module(VERIFY,'old_exact_fan_verifier_c012');verify.ROOT=HERE
        verify.RESULT=GEOMETRY;verify.main()
        replay=json.loads((HERE/'verify_exact_fan_results.json').read_text())
        assert replay==json.loads(OLD_VERIFY.read_text())
        old=json.loads(GEOMETRY.read_text());cones=parse_fan(OLD_FAN)
        assert len(cones)==len(old['cells'])==21
        V=[(1,)+s for s in product((-1,1),repeat=5)]
        axes=[tuple(int(i==j)for i in range(6))for j in range(6)]
        frame=V+axes
        A=fmpq_mat([[v[i]for v in frame]for i in range(6)])
        gram=A*A.transpose()
        assert all(gram[i,j]==(33 if i==j else 0)for i in range(6)for j in range(6))
        rows=[];text=['21'];total=Q(0);base_exp=Q(0);unique=set()
        for ci,((det,rays),cell)in enumerate(zip(cones,old['cells'])):
            matrix=fmpq_mat([[ray[i]for _,ray in rays]for i in range(6)])
            assert abs(matrix.det())==det
            scales=[];new_rates=[];text.append(str(det))
            for ri,((q,ray),original,qoriginal)in enumerate(zip(rays,cell['rays_a0_equals_one'],cell['rates'])):
                assert ray[0]>0 and all(ray[i]>=ray[i+1]>=0 for i in range(5))
                assert [Q(v,ray[0])for v in ray]==list(map(Q,original))
                assert Q(q,ray[0])==Q(qoriginal)
                body=sum(abs(sum(a*b for a,b in zip(v,ray)))for v in V)
                full=sum(abs(sum(a*b for a,b in zip(v,ray)))for v in frame)
                qnew=q+sum(map(abs,ray))
                assert q==body and qnew==full and qnew>q
                rows.append({'cone':ci,'ray_index':ri,'ray':ray,'q_old':q,
                             'axis_l1':sum(map(abs,ray)),'q_new':qnew,
                             'full_38_symbol_exact':full})
                unique.add(ray);scales.append(ray[0]);new_rates.append(qnew)
                text.append(' '.join(map(str,[qnew,*ray])))
            for v in frame:
                vals=[sum(a*b for a,b in zip(v,ray))for _,ray in rays]
                assert min(vals)>=0 or max(vals)<=0
            assert Q(det,math.prod(scales))==Q(cell['determinant'])
            total+=Q(det,math.prod(scales)*math.factorial(5))
            base_exp+=Q(det,math.prod(new_rates))
        assert total==Q(1,120)
        NEW_FAN.write_text('\n'.join(text)+'\n')
        native_record=json.loads((ENGINE/'build_record.json').read_text())
        native_check=json.loads((ENGINE/'validation_results_v3.json').read_text())
        assert sha(BINARY)==native_record['binary_sha256']==native_check['binary_sha256']
        assert sha(ENGINE/'arb_kernel.cpp')==native_check['native_source_sha256']
        dep=json.loads((ENGINE/'dependencies.json').read_text())
        for path,expected in dep['runtime_library_sha256'].items():assert sha(path)==expected
        # Input range safety is verified without converting or rounding either point.
        bounds=[]
        for point in points:
            ns=point['numerators'];den=point['denominator']
            assert ns[0]==den and all(ns[i]>ns[i+1]>0 for i in range(5))
            dot_bound=max(sum(abs(a)*max(map(abs,ns))for a in row['ray'])for row in rows)
            den_bound=den*max(row['q_new']for row in rows)
            safe=max(dot_bound,den_bound,max(map(abs,ns)),den)<2**63
            assert safe==(point['point_index']==1)
            bounds.append({'point_index':point['point_index'],'max_dot_numerator_bound':dot_bound,
                           'max_node_denominator_bound':den_bound,'native_int64_safe':safe,
                           'native_cross_product_int128_safe':dot_bound*den_bound<2**127})
        assert T/4096==Q(1,92160) and T/4000==Q(1,90000)
        assert Q(6,125)/T==Q(27,25)
        out={'status':'PASS_EXACT_AUGMENTED_FRAME_AND_FAN','pid':os.getpid(),
             'beta':'1','T':str(T),'columns':38,'Gram':'33 I6','old_cover_exact_replayed':True,
             'base_cones':21,'full_signed_permutation_group':46080,'complete_group_terms':967680,
             'ray_occurrences_checked':len(rows),'unique_integer_rays':len(unique),
             'section_volume':str(total),'all_coordinate_forms_nonnegative_on_each_base_cone':True,
             'all_38_linear_forms_constant_weak_sign_per_cone':True,'rays_determinants_unchanged':True,
             'rate_checks':rows,'new_base_exp_minus_symbol_integral':str(base_exp),
             'native_input_range_checks':bounds,'new_fan_sha256':sha(NEW_FAN),
             'source_sha256':sha(__file__),'manifest_sha256':sha(HERE/'input_manifest.json'),
             'fixed_T_note':{'explicit_T_kept':str(T),'old_reference_T':'6/125',
                 'actual_time_ratio':str(Q(6,125)/T),
                 '6_over_125_divided_by_1_plus_1_over_12':str(Q(6,125)/(1+Q(1,12))),
                 'parenthetical_scaling_is_not_an_exact_identity':True},
             'new_kernel_passes':0,'elapsed_seconds':perf_counter()-start}
        write('fan_audit.json',out)
        print(json.dumps({k:v for k,v in out.items()if k not in ['rate_checks']},indent=2))
    finally:signal.alarm(0)

if __name__=='__main__':main()
