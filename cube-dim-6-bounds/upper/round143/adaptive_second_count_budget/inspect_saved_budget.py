#!/usr/bin/env python3
"""C143 saved-data exact budget inspection; zero functional evaluation."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
import gzip,hashlib,json,sys,time
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
C130=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes')
C128=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round128/true_beta1_joint_moment')
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dec(q):
    with localcontext()as c:c.prec=24;return str(Decimal(q.numerator)/Decimal(q.denominator))
def encode(q):return {'exact':str(q),'display_only':dec(q)}
def main():
    start=time.monotonic();r=read(C130/'results.json');cones=read(C130/'cone_table.json')['cones'];spec=read(C130/'inputs.json')
    cert=read(C128/'candidate_certificate.json');candidate=read(C128/'candidate.json');upstream=read(C128/'upstream.json')
    lo,hi=map(F,r['moment_interval']);width=F(r['weighted_box_width']);tail=F(r['tail_upper']);piwidth=F(r['pi_factor_width'])
    assert hi-lo==width+tail+piwidth
    beta=F(spec['old_beta2_interval'][1]);assert beta==F(upstream['beta_by_h'][5][1])
    factorU=F(r['pi_factor_interval'][1]);leafrows=[];cone_stats=[];old_vertex_count=0;max_kernel_radius_width=F(0)
    sources=[C130/n for n in ['results.json','inputs.json','cone_table.json','box_core.py','run_prepared.py','parent_release.json','source_manifest.json']]+[C128/n for n in ['candidate_certificate.json','candidate.json','upstream.json','joint_core.py']]
    for cone,source in zip(cones,r['leaf_files']):
        path=C130/source['path'];assert sha(path)==source['sha256'];sources.append(path)
        with gzip.open(path,'rt')as f:data=json.load(f)
        assert data['cone_id']==cone['id'] and len(data['leaves'])==256
        seen={};sumwidth=F(0)
        for leaf in data['leaves']:
            L,U=map(F,leaf['integral_interval']);w=factorU*F(cone['alpha'])*(U-L);sumwidth+=w
            leafrows.append({'cone_id':cone['id'],'prefix':leaf['prefix'],'weighted_width':str(w),'integral_lower':leaf['integral_interval'][0],'integral_upper':leaf['integral_interval'][1],'lower_t':leaf['lower_t'],'upper_t':leaf['upper_t']})
            for side in ['lower','upper']:
                t=tuple(leaf[side+'_t']);v=leaf['kernel_at_'+side+'_t'];payload=json.dumps(v,separators=(',',':')).encode();digest=hashlib.sha256(payload).hexdigest()
                if t in seen:assert seen[t]==digest
                else:
                    seen[t]=digest
                    for pair in v:
                        a,b=map(F,pair);assert 0<=a<=b;max_kernel_radius_width=max(max_kernel_radius_width,b-a)
        old_vertex_count+=len(seen);cone_stats.append({'cone_id':cone['id'],'unique_saved_vertex_vectors':len(seen),'weighted_width':encode(sumwidth)})
    assert sum(F(x['weighted_width'])for x in leafrows)==width
    leafrows.sort(key=lambda x:(-F(x['weighted_width']),x['cone_id'],x['prefix']))
    gap=hi-beta;cum=F(0);necessary_k=None;prefixes={}
    for i,row in enumerate(leafrows,1):
        cum+=F(row['weighted_width'])
        if i in [15,30,60,120,240,480,960,1920,3840,5376]:prefixes[str(i)]=encode(cum)
        if necessary_k is None and cum>gap:necessary_k=i
    assert necessary_k is not None
    x=F(candidate['numerators'][32],int(candidate['denominator']));assert x>0
    old6=F(cert['M_power6_upper']);theta=F(cert['threshold_M_power6_exact'])
    required=beta*(1-(1-theta/old6)/x)
    fixed_ray_best_with_known_lower=(1-x+x*lo/beta)*old6
    # Uniform endpoint uncertainty bound for all stored boxes: pair_sum is
    # at most31*pU-Lipschitz in each of32 nonnegative kernel coordinates.
    pU=F(spec['p_interval'][1]);mass=sum(F(c['alpha'])for c in cones)*(1-F(1,2**24))**6
    kernel_uncertainty_upper=factorU*mass*1984*pU*max_kernel_radius_width
    out={'status':'SAVED_ONLY_EXACT_BUDGET_INSPECTION_NO_FUNCTIONAL_CALLS','source_130_primary_locked_science':True,
      'h':'1/128','mu2_interval':[encode(lo),encode(hi)],'old_beta2_upper':encode(beta),
      'width_decomposition_exact_verified':True,'weighted_box_width':encode(width),'tail_upper':encode(tail),'pi_factor_width':encode(piwidth),
      'gap_to_improve_old_beta2':encode(gap),'required_fraction_of_box_width_recovery':encode(gap/width),
      'top_weighted_width_prefix_sums':prefixes,'necessary_minimum_old_leaves_to_patch_even_with_ideal_width_recovery':necessary_k,
      'unique_recoverable_saved_vertex_vectors':old_vertex_count,'unique_recoverable_saved_scalar_kernel_records':32*old_vertex_count,
      'original_kernel_calls':r['kernel_calls'],'all_old_vertices_recoverable_proof':'Every split preserves both parent diagonal vertices in descendants; all evaluated diagonal vertices therefore occur among final leaf endpoints. Original root reset/LRU may have repeated calls.',
      'max_saved_kernel_interval_width':encode(max_kernel_radius_width),'global_endpoint_uncertainty_upper':encode(kernel_uncertainty_upper),
      'geometric_width_lower_after_endpoint_uncertainty':encode(max(F(0),width-kernel_uncertainty_upper)),
      'fixed128_ray':{'target_entry_zero_based':32,'simplex_weight_x':encode(x),'normalization_Z':'1-x+x*new_beta2/old_beta2','new_budget_power6':'Z*old_MU_power6, holding all unnormalized coefficients fixed',
        'beta2_threshold_for_g0':encode(required),'known_true_mu2_lower_exceeds_required':lo>required,
        'best_power6_with_known_mu2_lower_exceeds_g0':fixed_ray_best_with_known_lower>theta,
        'scope':'Only the current128 coefficient ray with the h1/128 beta2 replacement; no reoptimized-family exclusion'},
      'cone_stats':cone_stats,'zero_new_kernel_erfc_exp_log_or_root_calls':True,'old_functional_replays':0,
      'source_hashes':[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size}for p in sources],
      'script_sha256':sha(Path(__file__)),'wall_seconds':time.monotonic()-start}
    (P/'saved_budget_analysis.json').write_text(json.dumps(out,indent=2)+'\n')
    (P/'ranked_old_leaf_widths.json').write_text(json.dumps(leafrows,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ['status','gap_to_improve_old_beta2','required_fraction_of_box_width_recovery','necessary_minimum_old_leaves_to_patch_even_with_ideal_width_recovery','unique_recoverable_saved_vertex_vectors','unique_recoverable_saved_scalar_kernel_records','fixed128_ray','wall_seconds']},indent=2))
if __name__=='__main__':main()
