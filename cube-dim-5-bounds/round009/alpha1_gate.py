"""Exact necessary I_face gate; consume full source endpoints, never rerun kernels."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product, combinations
from math import isqrt
import argparse, gzip, hashlib, json, os, time

HERE = Path(__file__).resolve().parent
P = HERE.parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def root4(x, digits=100):
    den = 10**digits
    k = isqrt(isqrt((x.numerator*den**4)//x.denominator))
    lo, hi = F(k, den), F(k+1, den)
    assert lo**4 <= x < hi**4
    return {"lower_fraction": str(lo), "upper_fraction": str(hi)}
def decimal(x, digits=60):
    # Display only; all comparisons below use Fraction.
    scale = 10**digits
    n = x.numerator*scale//x.denominator
    return f"{n//scale}.{n%scale:0{digits}d}"
def in_pair_span(x, v, w):
    # Every distinct body pair has a nonzero 2x2 minor.
    i,j = next((i,j) for i,j in combinations(range(5),2)
               if v[i]*w[j]-v[j]*w[i])
    d = v[i]*w[j]-v[j]*w[i]
    a = (x[i]*w[j]-x[j]*w[i])/d
    b = (v[i]*x[j]-v[j]*x[i])/d
    return all(x[k] == a*v[k]+b*w[k] for k in range(5))
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-project', type=Path, default=P)
    parser.add_argument('--output', type=Path, default=HERE/'alpha1_gate.json')
    args = parser.parse_args(); S=args.source_project.resolve(); began=time.monotonic()
    green = S/'astra_green018/batch_result.json'
    table = S/'astra_vertex017/boundary_triangulation.json.gz'
    budget = P/'upper_alpha1_029/check_results.json'
    a=json.loads(green.read_text()); selected=a['records'][a['best_record']]
    assert selected == a['selected_value'] and a['status']=='complete'
    x=tuple(map(F,selected['point'])); assert x[0]==1 and all(0<z<1 for z in x[1:])
    gl,gu=F(selected['G_lower_dyadic']),F(selected['G_upper_dyadic'])
    assert 0<gl<=gu<F('122384/1000000000')
    t=json.loads(gzip.decompress(table.read_bytes())); assert t['denominator']==144
    assert len(t['vertices'])==280 and len(t['simplices'])==3504
    dots=[sum((F(z)*y for z,y in zip(v,x)),F(0))/144 for v in t['vertices']]
    assert all(dots) and len(set(dots))==280
    delta=min(map(abs,dots)); max_l1=max(sum(map(abs,v)) for v in t['vertices'])/F(144)
    radius=delta/(2*max_l1)  # infinity ball keeps every vertex projection nonzero
    V=[(1,)+s for s in product((-1,1),repeat=4)]
    pair_hits=[(i,j) for i,j in combinations(range(16),2) if in_pair_span(x,V[i],V[j])]
    assert pair_hits==[]
    b=json.loads(budget.read_text()); bounds=b['independent_512bit_checks']['I']['independent_interval']
    il,iu=F(bounds['lower_fraction']),F(bounds['upper_fraction'])
    cap=F(b['usable_strict_rational_I_upper']); assert 0<il<iu<cap
    necessary4=F(32**4,3125)*gu**5
    gate=root4(necessary4)
    cost4=3125*il**4/gu**5; cost=root4(cost4)
    residual=3125*il**4-32**4*gu**5; assert residual>0
    assert F(gate['upper_fraction'])<F('552/10000000')
    assert gu<F('147/1000000')  # also strictly improves saved old point obstacle
    sources=[Path(__file__),HERE/'alpha1_gate_registration.md',green,table,budget,
             P/'upper_alpha1_029/report.md',P/'audit_vertex032/report.md']
    out={'status':'complete_exact_same_model_budget_incompatibility_no_C5_bound',
         'model':'alpha1, 16 unnormalised body directions, same I_face as029',
         'point':[str(z) for z in x], 'G_point_full_dyadic_interval':[str(gl),str(gu)],
         'G_point_upper_display':decimal(gu), 'all_280_projections_nonzero_and_distinct':True,
         'minimum_abs_projection':str(delta), 'max_vertex_l1':str(max_l1),
         'nonzero_projection_neighbourhood_infinity_radius':str(radius),
         'body_pair_planes_checked':120,'body_pair_plane_hits':pair_hits,
         'essential_cube_g_upper':str(gu),
         'necessary_I_face_for_cost_below32_fourth_power':str(necessary4),
         'necessary_I_face_threshold_interval':gate,
         'necessary_I_face_threshold_upper_display':decimal(F(gate['upper_fraction'])),
         'old_complete_I_budget_interval':[str(il),str(iu)],'old_rational_cap':str(cap),
         'best_possible_cost_using_old_formula_lower_interval':cost,
         'best_possible_cost_lower_display':decimal(F(cost['lower_fraction'])),
         'strict_incompatibility_residual':str(residual),
         'normalization':'I_face=integral_[0,1]^4 G(1,y)^(5/4); positive-part mass=128 I_face',
         'limitations':['Not an upper bound on actual I_face','No global g lower bound',
                        'No C5 lower or upper improvement','Point is not claimed a global minimizer'],
         'source_sha256':{str(p):sha(p) for p in sources},
         'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began,'root_decimal_digits':100,
                    'production_kernel_evaluations':0,'new_probability_or_spatial_samples':0}}
    assert time.monotonic()-began<1
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','G_point_upper_display',
       'necessary_I_face_threshold_upper_display','best_possible_cost_lower_display','runtime']},indent=2))
if __name__=='__main__': main()
