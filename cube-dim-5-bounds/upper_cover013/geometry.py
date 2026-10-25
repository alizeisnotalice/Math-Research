#!/usr/bin/env python3
"""One exact finite geometry batch; no maximal-level-set sampling."""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import argparse
import hashlib
import json
import time

import sympy as sp


def volume(box):
    answer=F(1)
    for low,high in box: answer*=max(F(0),high-low)
    return answer
def intersect(a,b):
    return tuple((max(x[0],y[0]),min(x[1],y[1])) for x,y in zip(a,b))
def inside(point,box):
    return all(low<z<high for z,(low,high) in zip(point,box))
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--literature-dir',type=Path,help='Optional read-only frozen source directory; not a numerical dependency')
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args()
    started=time.monotonic();root=Path(__file__).resolve().parent
    registration=json.loads((root/'registration.json').read_text())
    R=sp.Symbol('R',positive=True)
    derivative=sp.factor(sp.diff((R-1)/R**5,R))
    assert sp.simplify(derivative-(5-4*R)/R**6)==0
    eta=F(registration['strict_witness_excess']);d=F(registration['receiver_neighborhood_radius'])
    eta_core=F(registration['core_mass_excess_for_union_obstruction'])
    assert (1+eta)/(1+d)**5>1
    split=F(2,9);core_constant=F(17,10)**5;local_bound=F(256,9)
    assert F(6,5)**5>2 and core_constant<F(128,9)
    assert 32-16*split==local_bound and core_constant+64*split<local_bound<32
    signs=list(product((-1,1),repeat=5));V=F(32);rows=[]
    for text in registration['preselected_offsets']:
        a=F(text);assert F(1,2)<a<1
        centers=[tuple(a*s for s in sig) for sig in signs]
        cubes=[tuple((c-1,c+1) for c in center) for center in centers]
        central=tuple((F(-1),F(1)) for _ in range(5))
        core=tuple((a-1,1-a) for _ in range(5))
        pair_distance=min(max(abs(x-y) for x,y in zip(c1,c2)) for c1,c2 in combinations(centers,2))
        assert pair_distance==2*a>1 and all(inside((F(0),)*5,Q) for Q in cubes)
        endpoints=sorted({-1-a,F(-1),a-1,1-a,F(1),1+a})
        intervals=list(zip(endpoints,endpoints[1:]))
        assert len(intervals)==5
        in_hist={};out_hist={};N_integral=F(0);union=F(0);cells=0
        for cell in product(intervals,repeat=5):
            cells+=1;mid=tuple((x+y)/2 for x,y in cell);v=volume(cell)
            N=sum(inside(mid,Q) for Q in cubes)
            assert N>=1
            hist=in_hist if inside(mid,central) else out_hist
            hist[N]=hist.get(N,F(0))+v
            N_integral+=N*v;union+=v
        assert cells==3125 and N_integral==32*V
        assert union==(2*(1+a))**5
        assert max(in_hist)==32 and in_hist[32]==volume(core)
        assert max(out_hist)==16 and out_hist[16]>0
        # Ten disjoint outer face boxes; every Q contains exactly one per axis.
        slabs=[]
        for j in range(5):
            for sign in (-1,1):
                slab=list(core)
                slab[j]=(F(1),1+a) if sign==1 else (-1-a,F(-1))
                slabs.append(tuple(slab))
        slab_volume=a*(2*(1-a))**4
        assert all(volume(S)==slab_volume and volume(intersect(S,central))==0 for S in slabs)
        assert all(volume(intersect(S,T))==0 for S,T in combinations(slabs,2))
        incidence=[[volume(intersect(Q,S))/slab_volume for S in slabs] for Q in cubes]
        assert all(set(row)<={F(0),F(1)} and sum(row)==5 for row in incidence)
        mass_per_slab=(1+eta)*V/5;total_mass=10*mass_per_slab
        cube_integrals=[mass_per_slab*sum(row) for row in incidence]
        assert all(m/V==1+eta for m in cube_integrals)
        assert all(sum(row[j] for row in incidence)==16 for j in range(10))
        kernel_integral_against_f=sum(cube_integrals)/32
        assert kernel_integral_against_f==total_mass/2==(1+eta)*V
        peak_radius=F(5,4)
        assert 1<peak_radius<1+a
        peak=(1+eta)*F(512,3125)/a
        assert peak==2*(1+eta)*(peak_radius-1)/(a*peak_radius**5)<1
        assert pair_distance>2*d
        neighborhood_volume=32*(2*d)**5
        core_input_mass=(1+eta_core)*V
        union_mass_ratio=union/core_input_mass
        core_receiver_halo_ratio=(2-a+eta_core)**5/(1+eta_core)
        if a==F(999,1000):
            assert union_mass_ratio>F(319,10) and core_receiver_halo_ratio<F(101,100)
        rows.append({
            'a':str(a),'cube_count':32,'radius':'1','pair_center_linf_distance':str(pair_distance),
            'cell_count':cells,'common_core_volume':str(volume(core)),
            'union_volume':str(union),'sum_cube_volumes':str(N_integral),
            'union_over_sum_volume':str(union/N_integral),
            'inside_central_multiplicity_volume':{str(k):str(v) for k,v in sorted(in_hist.items())},
            'outside_central_multiplicity_volume':{str(k):str(v) for k,v in sorted(out_hist.items())},
            'outside_multiplicity_max':16,'packet_kernel_integral':'32',
            'outer_input':{'slab_count':10,'slab_volume':str(slab_volume),
                'constant_height_on_each_slab':str(mass_per_slab/slab_volume),
                'total_mass':str(total_mass),'mass_in_central_cube':'0',
                'every_witness_cube_average':str(1+eta),
                'common_center_maximal_exact':str(peak),'common_center_maximizing_radius':'5/4',
                'common_center_not_in_strict_level_set':True,
                'kernel_integral_against_f':str(kernel_integral_against_f),
                'receiver_neighborhood_radius':str(d),
                'receiver_neighborhood_union_volume':str(neighborhood_volume),
                'receiver_neighborhood_average_lower':str((1+eta)/(1+d)**5)},
            'core_supported_input':{'mass':str(core_input_mass),
                'every_witness_average':str(1+eta_core),
                'selected_cube_union_volume_over_input_mass':str(union_mass_ratio),
                'actual_level_set_volume_over_mass_strict_upper':str(core_receiver_halo_ratio),
                'warning':'Selected union is not asserted to lie in the maximal strict level set.'}})
    sources={p.name:sha(p) for p in [Path(__file__),root/'registration.json']}
    if args.literature_dir:
        for name in ['strict_upper_second_batch.md','table.md']:
            sources['frozen_literature/'+name]=sha(args.literature_dir/name)
    payload={'status':'finite_local_lemmas_and_counterexamples_verified_not_global_upper_improvement',
        'C5_improved':False,'source_URL':'https://arxiv.org/html/1811.10478v2',
        'source_read':'Definition1.1, Theorem1.2/2.5, Lemma2.6 proof chain, Theorem3.3, Corollary3.5',
        'symbolic_derivative':str(derivative),'rows':rows,'source_sha256':sources,
        'uniform_local_star_lemma':{'outside_core_mass_fraction_split':str(split),
            'core_receiver_bound':str(core_constant),'outside_weak_type_charge_factor':'64',
            'high_outside_mass_bound_at_split':str(32-16*split),
            'low_outside_mass_bound_at_split':str(core_constant+64*split),
            'universal_local_bound':str(local_bound),
            'saving_from32':str(32-local_bound),
            'receiver':'{Mf>lambda} intersect union of the 32 qualifying star cubes',
            'RHS_mass':'full ||f||_1; not additive across arbitrary packets'},
        'runtime':{'elapsed_seconds':time.monotonic()-started,'sympy_version':sp.__version__,
                   'arithmetic':'Fraction exact rational geometry; one SymPy derivative identity',
                   'cells_examined':3*3125,'maximal_function_sample_grid_points':0},
        'remaining_gap':'No global packet selection covering every compact strict level set with pointwise weighted congestion <32 has been constructed.'}
    args.output_dir.mkdir(parents=True,exist_ok=True)
    destination=args.output_dir/'geometry_results.json'
    destination.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({'output':str(destination),'sha256':sha(destination),
                      'elapsed_seconds':payload['runtime']['elapsed_seconds'],
                      'examples':[{'a':r['a'],
                        'common_center_maximal_exact':r['outer_input']['common_center_maximal_exact'],
                        'all_witness_averages':r['outer_input']['every_witness_cube_average'],
                        'union_core_mass_ratio':r['core_supported_input']['selected_cube_union_volume_over_input_mass']}
                        for r in rows]},indent=2))


if __name__=='__main__': main()
