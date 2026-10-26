#!/usr/bin/env python3
"""One small Fraction qualification, no union and no baseline evaluation."""
import time
START=time.monotonic()
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib,json,os,sys
import generator as gen

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    reg=json.loads((HERE/'registration.json').read_text())
    assert reg['status']=='registered_before_small_fraction_qualification'
    assert reg['source_sha256']==sha(__file__)
    assert reg['generator_sha256']==sha(HERE/'generator.py')
    assert not (HERE/'results.json').exists()
    out={'pid':os.getpid(),'python':sys.executable,'status':'running',
         'registration_sha256':sha(HERE/'registration.json'),'source_sha256':sha(__file__),
         'generator_sha256':sha(HERE/'generator.py')}
    print(json.dumps(out),flush=True)
    try:
        for item in reg['readonly_inputs']:
            assert sha(item['path'])==item['sha256'],item['path']
        assert Q(9,8)**5<2<Q(6,5)**5
        assert Q(3,2)**5<8<Q(8,5)**5
        assert Q(65,64)**5>Q(17,16)
        rmax=Q(65,32);rsingle=Q(13,8)
        assert Q(15,16)*Q(7,4)>rsingle
        assert Q(31,16)*Q(7,4)>rmax
        assert Q(31,16)*Q(17,16)>rmax
        residuals={'first_singleton':Q(15,16)*Q(7,4)-rsingle,
                   'first_low_or_multicolumn_high':Q(31,16)*Q(7,4)-rmax,
                   'second_both_sides':Q(31,16)*Q(17,16)-rmax}
        # Boundary-rational test spacings, not a new lattice candidate or C evaluation.
        a,s=Q(7,4),Q(17,16)
        cases=0;atoms_checked=0
        for j1,j2,l1,l2,se,sh,sm in product(range(2),range(2),range(1,5),range(1,5),(-1,0,1),(-1,0,1),(-1,0,1)):
            eps,eta,mu=se*a/16,sh*s/16,Q(sm,16)
            lengths=(l1,l2,1,1,1)
            indices=[(p,q,0,0,0) for p in range(j1,j1+l1) for q in range(j2,j2+l2)]
            actual=[gen.atom(j,a,s,eps,eta,mu) for j in indices]
            assert len(set(x for x,w in actual))==l1*l2
            assert all(w>0 for x,w in actual)
            claimed=gen.extrema(lengths,(j1,j2),a,s,eps,eta)
            observed=tuple((min(x[k] for x,w in actual),max(x[k] for x,w in actual)) for k in range(5))
            assert claimed==observed
            assert sum(w for x,w in actual)==gen.mass(lengths,j1,mu)[0]
            for j,(x,w) in zip(indices,actual):
                reflected=(-j[0],j[1]+1,j[2],j[3],j[4])
                y,wy=gen.atom(reflected,a,s,eps,eta,mu)
                assert y==(-x[0],x[1]+s,*x[2:]) and wy==w
            cases+=1;atoms_checked+=len(actual)
        assert cases==1728
        assert sum(1 for _ in gen.descriptors())==16384
        # The e/h/m endpoint interface is checked algebraically at a formal r0=1.
        # This is a derivative-formula check, not an evaluation of nominal N^(1/5)/2.
        jet_cases=0
        for l1,l2,j1,j2 in product(range(1,5),range(1,5),range(2),range(2)):
            lengths=(l1,l2,1,1,1)
            row=gen.raw_receiver_jet(lengths,(j1,j2),a,s,Q(1),(a,s,Q(1)))
            expected=Q(gen.parity(j1),5*l1) if l1%2 else Q(0)
            assert row['radius_jet']==(Q(1),expected)
            assert row['mass_speed']==(gen.parity(j1)*l2 if l1%2 else 0)
            for k in (2,3,4):
                assert row['extrema_speed'][k]==(0,0)
            jet_cases+=1
        # Explicit period-index collar after the upper-half glide and tail reflections.
        q1=set();q2=set()
        for j in range(-1,6):
            q1.update((j//2,(2-j)//2));q2.update((j//2,(j+1)//2))
        assert (min(q1),max(q1))==(-2,2)
        assert (min(q2),max(q2))==(-1,3)
        for ell in range(1,5):
            indices=gen.tail_indices(ell)
            assert min(indices)>=-1 and max(indices)<=2
            assert all(-1<=1-j<=2 for j in indices)
        for item in reg['readonly_inputs']:
            assert sha(item['path'])==item['sha256'],item['path']
        out.update(status='qualified_geometry_and_interface_only_no_union',active_cases=cases,
                   actual_distinct_atom_rows_checked=atoms_checked,jet_formula_cases=jet_cases,
                   complete_raw_descriptor_count=16384,
                   strict_rational_start_guard_residuals={k:str(v) for k,v in residuals.items()},
                   collar_period_ranges={'first':[-2,2],'second':[-1,3],'tails':[-1,2],
                                         'uniform_width':3},
                   period_mass=4,full_period_receiver_copy_factor=16,ratio_multiplier=4,
                   new_union_calls=0,old_baseline_calls=0,new_parameter_grid=False,
                   optimality_or_improvement_claim=False)
    except Exception as error:
        out.update(status='stopped_failure_no_rerun',error=repr(error))
        raise
    finally:
        out['seconds_before_final_write']=time.monotonic()-START
        (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','pid','active_cases','actual_distinct_atom_rows_checked','seconds_before_final_write')}),flush=True)


if __name__=='__main__':main()
