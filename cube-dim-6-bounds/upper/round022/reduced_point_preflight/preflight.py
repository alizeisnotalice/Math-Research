#!/usr/bin/env python3
"""Exact one-point orbit/tuple compression only; no kernel evaluation."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter, defaultdict
from itertools import permutations, product
from time import perf_counter
import hashlib,json,math,os,resource,signal,sys

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
FAN=PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
GEOMETRY=PROJECT/'upper/round002/fan/exact_cells_results.json'
EXPECTED_FAN='7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'
X=(50,50,50,45,27,0);DEN=50
MAX_SECONDS=1200;MAX_RSS=768*1024**2;CEILING=1024**3

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(name,value):(HERE/name).write_text(json.dumps(value,indent=2)+'\n')
def peak_rss():
    n=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(n if sys.platform=='darwin'else n*1024)
def negate(nodes):return tuple(-v for v in reversed(nodes))
def parse_fan():
    tokens=iter(FAN.read_text().split());cones=[]
    for _ in range(int(next(tokens))):
        det=int(next(tokens));rays=[];alpha=Q(det)
        for _ in range(6):
            rate=int(next(tokens));ray=tuple(int(next(tokens))for _ in range(6))
            assert rate>0;rays.append((rate,ray));alpha/=rate
        cones.append((alpha,rays))
    assert next(tokens,None)is None and len(cones)==21
    return cones

def main():
    start=perf_counter();logical=0;actual=0
    sources=[FAN,GEOMETRY,PROJECT/'upper/round002/fan/verify_exact_fan.py',
             PROJECT/'upper/round002/fan/verify_exact_fan_results.json',
             PROJECT/'upper/round017/time_family_obstruction/theorem.md',
             HERE/'registration.md',HERE/'theorem.md',Path(__file__)]
    assert sha(FAN)==EXPECTED_FAN
    manifest={'status':'REGISTERED_BEFORE_STRUCTURE_ENUMERATION',
              'point_numerators':X,'denominator':DEN,'max_logical_terms':967680,
              'max_seconds':MAX_SECONDS,'memory_ceiling_bytes':CEILING,'rss_guard_bytes':MAX_RSS,
              'inputs':[{'path':str(p),'sha256':sha(p)}for p in sources]}
    write('input_manifest.json',manifest)
    def guard():
        if perf_counter()-start>MAX_SECONDS:raise TimeoutError('registered 1200-second cap')
        if peak_rss()>=MAX_RSS:raise MemoryError('conservative 768-MiB RSS guard')
        if logical>967680:raise RuntimeError('logical term cap')
    signal.signal(signal.SIGALRM,lambda a,b:(_ for _ in ()).throw(TimeoutError('registered 1200-second cap')))
    signal.alarm(MAX_SECONDS)
    try:
        orbit=Counter();stabilizer=[];group_count=0
        for perm in permutations(range(6)):
            for mask in range(64):
                y=tuple((-1 if mask>>j&1 else 1)*X[perm[j]]for j in range(6))
                orbit[y]+=1;group_count+=1
                if y==X:stabilizer.append({'permutation':perm,'sign_mask':mask})
            guard()
        assert group_count==46080 and len(stabilizer)==12 and len(orbit)==3840
        assert set(orbit.values())=={12} and sum(orbit.values())==46080
        cones=parse_fan();groups={};node_values=set();squared_values=set()
        zero_terms=repeat_terms=both_terms=0;zero_hist=Counter();multiplicity_hist=Counter()
        max_factor_mult=defaultdict(int);total_weight=Q(0)
        for oi,(y,multiplicity)in enumerate(orbit.items()):
            for alpha,rays in cones:
                nodes=tuple(sorted(Q(sum(y[j]*ray[j]for j in range(6)),DEN*q)for q,ray in rays))
                logical+=multiplicity;actual+=1
                counts=Counter(nodes);zero_count=counts.get(Q(0),0)
                repeated=len(counts)<6
                zero_terms+=multiplicity*bool(zero_count)
                repeat_terms+=multiplicity*repeated
                both_terms+=multiplicity*(bool(zero_count)and repeated)
                zero_hist[zero_count]+=multiplicity
                multiplicity_hist[tuple(sorted(counts.values(),reverse=True))]+=multiplicity
                squares=Counter(v*v for v in nodes)
                for value,power in squares.items():max_factor_mult[value]=max(max_factor_mult[value],power)
                node_values.update(nodes);squared_values.update(squares)
                weight=alpha*multiplicity;total_weight+=weight
                entry=groups.setdefault(nodes,[Q(0),0,0])
                entry[0]+=weight;entry[1]+=multiplicity;entry[2]+=1
            if oi%32==0:guard()
        assert logical==967680 and actual==80640
        assert total_weight==Q(json.loads(GEOMETRY.read_text())['full_exp_minus_symbol_integral'])
        assert sum(v[1]for v in groups.values())==logical
        assert sum((v[0]for v in groups.values()),Q(0))==total_weight
        canonical={};self_pairs=0
        for nodes,(weight,count,orbit_terms)in groups.items():
            other=negate(nodes)
            assert other in groups and groups[other][0]==weight and groups[other][1]==count
            key=min(nodes,other)
            if nodes==other:self_pairs+=1
            entry=canonical.setdefault(key,[Q(0),0,0])
            entry[0]+=weight;entry[1]+=count;entry[2]+=orbit_terms
        assert sum(v[1]for v in canonical.values())==logical
        assert sum((v[0]for v in canonical.values()),Q(0))==total_weight
        assert len(groups)==2*len(canonical)-self_pairs
        guard()
        # Exact geometry check from four selected coordinate forms, not frame enumeration.
        selected=[Q(X[i],DEN)for i in (0,3,4,5)]
        assert len(set(map(abs,selected)))==4
        pair_values=[abs(selected[i]+s*selected[j])for i in range(4)for j in range(i)for s in (-1,1)]
        signed_values=[abs(selected[0]+sum(s*v for s,v in zip(signs,selected[1:])))for signs in product((-1,1),repeat=3)]
        assert min(pair_values)==Q(1,10)and min(signed_values)==Q(11,25)
        distance_lower_sq=min(min(pair_values)**2/2,min(signed_values)**2/4)
        assert distance_lower_sq==Q(1,200)
        merged=[]
        for nodes,(weight,count,orbit_terms)in sorted(canonical.items()):
            merged.append({'nodes':list(map(str,nodes)),'weight':str(weight),
                           'logical_terms':count,'point_orbit_cone_terms':orbit_terms,
                           'global_negation_orbit_size':1 if nodes==negate(nodes)else 2})
        write('merged_tuples.json',{'status':'EXACT_FULL_WEIGHT_NODE_TUPLE_COMPRESSION_ONLY',
             'point_numerators':X,'denominator':DEN,'logical_terms':logical,
             'tuples':merged,'density_evaluated':False})
        factors=[{'b_squared':str(value),'raw_exponent':power}for value,power in sorted(max_factor_mult.items())]
        raw_degree=sum(max_factor_mult.values())
        write('denominator_factors.json',{'variable':'z=T^2','factor_form':'z+b_squared',
             'factors':factors,'raw_degree_z':raw_degree,'raw_degree_T':2*raw_degree,
             'no_numerator_formed':True,'no_cancellation_performed':True,'reduced_degree':None})
        guard()
        out={'status':'PASS_COMPLETE_STRUCTURAL_PREFLIGHT_NO_KERNEL',
             'pid':os.getpid(),'python':sys.version,'point_numerators':X,'denominator':DEN,
             'group_elements_enumerated':group_count,'stabilizer_size':len(stabilizer),
             'stabilizer_elements':stabilizer,'distinct_point_orbit_images':len(orbit),
             'uniform_image_multiplicity':12,'actual_orbit_cone_pairs':actual,
             'logical_full_group_terms':logical,'distinct_exact_b':len(node_values),
             'distinct_exact_b_squared':len(squared_values),
             'minimum_abs_nonzero_b':str(min(abs(v)for v in node_values if v)),
             'zero_node_logical_terms':zero_terms,'repeated_node_logical_terms':repeat_terms,
             'both_zero_and_repeated_logical_terms':both_terms,
             'zero_multiplicity_histogram':{str(k):v for k,v in sorted(zero_hist.items())},
             'node_multiplicity_pattern_histogram':{' '.join(map(str,k)):v for k,v in sorted(multiplicity_hist.items())},
             'unordered_signed_node_tuples':len(groups),'tuples_mod_global_negation':len(canonical),
             'self_negation_symmetric_tuples':self_pairs,
             'squared_node_multiset_count':len({tuple(sorted(v*v for v in nodes))for nodes in groups}),
             'total_full_group_alpha_weight':str(total_weight),
             'max_raw_pole_multiplicity':max(max_factor_mult.values()),
             'raw_factor_exponent_histogram':dict(sorted(Counter(max_factor_mult.values()).items())),
             'zero_pole_raw_exponent':max_factor_mult.get(Q(0),0),
             'raw_denominator_degree_z':raw_degree,'raw_denominator_degree_T':2*raw_degree,
             'raw_upper_bound_within_prior_512_gate':raw_degree<=512,
             'reduced_denominator_degree':None,'numerator_formed':False,'cancellation_performed':False,
             'off_S2_proved':True,'off_S3_proved':True,
             'selected_coordinates':list(map(str,selected)),
             'selected_min_pair_form_abs':str(min(pair_values)),
             'selected_signed_sum_abs_values':list(map(str,sorted(set(signed_values)))),
             'S3_distance_squared_certified_lower':str(distance_lower_sq),
             'generic_zero_excluding_engine_usable':zero_terms==0,
             'new_G_H_K_p_evaluations':0,'new_time_evaluations':0,'Sturm_started':False,
             'new_spatial_points_evaluated':0,'structural_points_processed':1,
             'all_T_claim':False,'source_sha256':sha(__file__),
             'input_manifest_sha256':sha(HERE/'input_manifest.json'),
             'merged_tuples_sha256':sha(HERE/'merged_tuples.json'),
             'denominator_factors_sha256':sha(HERE/'denominator_factors.json'),
             'peak_rss_bytes':peak_rss(),'elapsed_seconds':perf_counter()-start,'exit_code':0}
        assert out['peak_rss_bytes']<CEILING
        write('results.json',out)
        print(json.dumps({k:v for k,v in out.items()if k not in ['stabilizer_elements']},indent=2))
    except BaseException as ex:
        write('stopped.json',{'status':'STOPPED_NO_FURTHER_WORK','error':repr(ex),
              'logical_terms_completed':logical,'actual_pairs_completed':actual,
              'peak_rss_bytes':peak_rss(),'elapsed_seconds':perf_counter()-start,
              'new_G_H_K_p_evaluations':0,'exit_code':1})
        raise
    finally:signal.alarm(0)

if __name__=='__main__':main()
