#!/usr/bin/env python3
"""One fmpq full-symbolic pass, no numerical time and no common numerator."""
import sys
sys.dont_write_bytecode = True
import argparse
import gzip
import hashlib
import json
import math
import os
import resource
import time
from collections import Counter, defaultdict
from pathlib import Path
from flint import fmpq as Q
import flint

HERE = Path(__file__).resolve().parent
SOURCE = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round022/reduced_point_preflight')
ARCHIVE = SOURCE/'merged_tuples.json.gz'
ARCHIVE_SHA = '562e3703286fc9a9712efdad30bbebb71433c4f3f15f867c985ee41b80c6f719'
CONTENT_SHA = '7deb058449e11fa74d7e1a4658311cb738a38fbbb57ab03118b494138cd39390'
REPORT_SHA = 'fad56fada09c08e96311fa7bf25f8d5326ca67947415688febb964d47bac72c5'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False)+'\n')

def read_input():
    compressed = ARCHIVE.read_bytes()
    assert hashlib.sha256(compressed).hexdigest() == ARCHIVE_SHA
    content = gzip.decompress(compressed)
    assert hashlib.sha256(content).hexdigest() == CONTENT_SHA
    data = json.loads(content)
    assert data['point_numerators'] == [50,50,50,45,27,0] and data['denominator'] == 50
    assert data['logical_terms'] == 967680 and not data['density_evaluated']
    assert len(data['tuples']) == 40320
    return data

def pole_coefficients(nodes):
    """Return real a_(b,j); complex coefficient is i^(j-6)*a_(b,j)."""
    counts = Counter(nodes)
    output = {}
    for b,m in counts.items():
        coefficients = [Q(1)] + [Q(0)]*(m-1)
        for c,n in counts.items():
            if c == b:
                continue
            delta = b-c
            factor = [Q((-1)**k*math.comb(n+k-1,k))*delta**(-n-k) for k in range(m)]
            coefficients = [sum((coefficients[h]*factor[k-h] for h in range(k+1)),Q(0)) for k in range(m)]
        for k,a in enumerate(coefficients):
            output[(b,m-k)] = a
    return output

def fold_positive(table):
    positive = sorted({abs(b) for b,j in table if b})
    folded = {}
    for b in positive:
        for j in range(1,7):
            c = table.get((b,j),Q(0))+(-1)**j*table.get((-b,j),Q(0))
            if c:
                folded[(b,j)] = c
    return folded

def numerator_coefficients(b,j):
    # Ascending powers of z; no coefficient depends on a sampled T.
    return [Q(1 if (j-r-3)%2 == 0 else -1)*math.comb(j,2*r)*b**(j-2*r)
            for r in range(j//2+1)]

def real_principal(folded):
    result = defaultdict(lambda:Q(0))
    for (b,j),c in folded.items():
        q = b*b
        polynomial = numerator_coefficients(b,j)
        for s in range(len(polynomial)):
            shifted = sum((polynomial[r]*math.comb(r,s)*(-q)**(r-s)
                           for r in range(s,len(polynomial))),Q(0))
            result[(q,j-s)] += c*shifted
    return {key:value for key,value in result.items() if value}

def rows(table, keyname):
    return [{keyname:str(b),'order':j,'coefficient':str(a)}
            for (b,j),a in sorted(table.items()) if a]

def register():
    assert not (HERE/'input_manifest.json').exists(), 'registration already exists'
    assert sha(SOURCE/'report.md') == REPORT_SHA
    data = read_input()
    files = [ARCHIVE,SOURCE/'report.md',SOURCE/'theorem.md',SOURCE/'results.json',
             SOURCE/'denominator_factors.json',SOURCE/'archive_receipt.json',
             HERE/'registration.md',HERE/'compress.py',HERE/'check.py',HERE/'theorem.md']
    manifest = {'status':'REGISTERED_BEFORE_SYMBOLIC_FULL', 'registered_epoch':time.time(),
        'pid':os.getpid(),'point_numerators':data['point_numerators'],'denominator':50,
        'symbolic_full_cap':1,'numeric_kernel_cap':0,'positive_time_values_cap':0,
        'tuple_cap':40320,'logical_terms_cap':967680,'scalar_case_cap':8,
        'seconds_cap':1200,'rss_bytes_cap':768*1024**2,'decompressed_sha256':CONTENT_SHA,
        'inputs':[{'path':str(p),'sha256':sha(p)}for p in files],
        'dependencies':{'python':sys.version,'python_flint':flint.__version__}}
    dump(HERE/'input_manifest.json',manifest)
    print(json.dumps({'status':manifest['status'],'pid':os.getpid(),'input_manifest_sha256':sha(HERE/'input_manifest.json')},indent=2))

def run():
    assert not (HERE/'run_receipt.json').exists(), 'one symbolic full pass only'
    manifest = json.loads((HERE/'input_manifest.json').read_text())
    for row in manifest['inputs']:
        assert sha(row['path']) == row['sha256']
    anchors = json.loads((HERE/'scalar_checks.json').read_text())
    assert anchors['status'] == 'PASS_EIGHT_FORMAL_PARTIAL_IDENTITIES'
    assert anchors['cases_count'] == 8 and anchors['numeric_kernel_calls'] == 0
    data = read_input()
    raw_factors = json.loads((SOURCE/'denominator_factors.json').read_text())['factors']
    raw_map = {Q(r['b_squared']):r['raw_exponent'] for r in raw_factors}
    assert len(raw_map) == 792 and sum(raw_map.values()) == 889
    started = time.perf_counter()
    table = defaultdict(lambda:Q(0))
    record = {'status':'ONE_SYMBOLIC_FULL_RUNNING','pid':os.getpid(),'tuples_processed':0,
              'logical_terms':0,'point_orbit_cone_terms':0,'symbolic_full_passes':1,
              'numeric_kernel_calls':0,'positive_time_values':0,'raw_partial_contributions':0,
              'zero_tuple_logical_terms':0,'repeated_tuple_logical_terms':0,
              'input_manifest_sha256':sha(HERE/'input_manifest.json'),'scalar_checks_sha256':sha(HERE/'scalar_checks.json')}
    dump(HERE/'run_receipt.json',record)
    weights = Q(0)
    try:
        for index,row in enumerate(data['tuples']):
            if index % 128 == 0:
                elapsed = time.time()-manifest['registered_epoch']
                rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                if sys.platform != 'darwin':
                    rss *= 1024
                assert elapsed < manifest['seconds_cap'], 'registered time cap'
                assert rss < manifest['rss_bytes_cap'], 'resident-memory cap'
                record['peak_rss_bytes'] = rss
            nodes = tuple(Q(s) for s in row['nodes'])
            weight = Q(row['weight'])
            assert len(nodes) == 6 and nodes == tuple(sorted(nodes)) and weight > 0
            assert row['global_negation_orbit_size'] == 2
            assert row['logical_terms'] == 12*row['point_orbit_cone_terms']
            counts = Counter(nodes)
            assert max(counts.values()) <= 4
            partial = pole_coefficients(nodes)
            for key,a in partial.items():
                table[key] += weight*a  # the stored weight already includes all multiplicities
            weights += weight
            record['tuples_processed'] += 1
            record['logical_terms'] += row['logical_terms']
            record['point_orbit_cone_terms'] += row['point_orbit_cone_terms']
            record['raw_partial_contributions'] += len(partial)
            record['zero_tuple_logical_terms'] += row['logical_terms']*int(Q(0) in counts)
            record['repeated_tuple_logical_terms'] += row['logical_terms']*int(max(counts.values())>1)
            if (index+1) % 2048 == 0:
                record['elapsed_seconds'] = time.perf_counter()-started
                dump(HERE/'run_receipt.json',record)
        assert record['tuples_processed'] == 40320 and record['logical_terms'] == 967680
        assert record['point_orbit_cone_terms'] == 80640
        assert weights == Q(520066141,48769479475200)
        assert record['zero_tuple_logical_terms'] == 88320 and record['repeated_tuple_logical_terms'] == 201408
        zero = [{'order':j,'A':str(table.get((Q(0),j),Q(0))),
                 'real_effect':'must cancel exactly' if j%2 == 0 else 'identically zero: i^(j-6) is imaginary'}for j in range(1,7)]
        assert all(table.get((Q(0),j),Q(0)) == 0 for j in (2,4,6)), 'uncancelled even zero pole'
        nonzero = {key:a for key,a in table.items() if a}
        folded = fold_positive(nonzero)
        principal = real_principal(folded)
        exponents = {}
        for b,j in folded:
            exponents[b*b] = max(exponents.get(b*b,0),j)
        actual_exponents = {}
        for q,j in principal:
            actual_exponents[q] = max(actual_exponents.get(q,0),j)
        assert exponents == actual_exponents
        assert all(q in raw_map and j <= raw_map[q] for q,j in exponents.items())
        gamma = sum((c for (b,j),c in folded.items() if j == 1),Q(0))
        assert gamma > 0
        output = {'status':'EXACT_COMPLETE_PARTIAL_TABLES','normalization':'p_T=(2*pi)^-6 R(T^2)',
                  'signed_A':rows(nonzero,'b'),'positive_fold_C':rows(folded,'positive_b'),
                  'real_z_principal':rows(principal,'b_squared'), 'zero_pole_channels':zero,
                  'denominator_factors':[{'b_squared':str(q),'exponent':j,'raw_exponent':raw_map[q]}
                                          for q,j in sorted(exponents.items())]}
        dump(HERE/'partial_tables.json',output)
        summary = {'status':'FULL_SYMBOLIC_COMPRESSION_COMPLETE',
            'signed_slots_touched':len(table),'nonzero_signed_slots':len(nonzero),
            'nonzero_positive_fold_slots':len(folded),'nonzero_real_principal_slots':len(principal),
            'raw_factor_count':792,'raw_denominator_degree':889,
            'reduced_factor_count':len(exponents),'reduced_denominator_degree':sum(exponents.values()),
            'reduced_exponent_histogram':dict(sorted(Counter(exponents.values()).items())),
            'zero_even_channels_exactly_cancelled':True,'zero_channels':zero,
            'full_weight_checksum':str(weights),'Gamma_G':str(gamma),'G_times_pi5_exact':str(gamma/128),
            'G_formula':'G(x0)=Gamma_G/(128*pi^5)',
            'symbolic_full_passes':1,'numeric_kernel_calls':0,'positive_time_values':0,
            'no_global_numerator_formed':True,'no_Sturm_or_root_isolation':True,
            'partial_tables_sha256':sha(HERE/'partial_tables.json')}
        dump(HERE/'results.json',summary)
        record['status'] = 'SYMBOLIC_FULL_COMPLETE_STAGE1_STOPPED'
        record['summary_counts'] = {k:v for k,v in summary.items() if 'count' in k or 'slots' in k or 'degree' in k}
    except BaseException as error:
        record['status'] = 'FAILED_STOPPED_NOT_A_CERTIFICATE'
        record['exception'] = repr(error)
        dump(HERE/'incomplete_partial_state.json',{'complete':False,'partial_signed_A':rows(table,'b')})
        raise
    finally:
        record['elapsed_seconds'] = time.perf_counter()-started
        dump(HERE/'run_receipt.json',record)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('Gamma_G','G_times_pi5_exact')},indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode',choices=('register','run'))
    args = parser.parse_args()
    (register if args.mode == 'register' else run)()
