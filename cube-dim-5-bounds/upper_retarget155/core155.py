"""Future153 joint target115 core. No real q/log/corner call on import."""
from fractions import Fraction as F
from itertools import product
import importlib.util
from util155 import sha
PENALTY_PATH='/Users/zhengzhihao/Desktop/T/Math-Research-cube5-curvature128-sol/cube-dim-5-bounds/upper_curvature153/penalty.py'
PENALTY_SHA='0f271bf3ffc6a52bb276d16a1144f4b0e748653054ce4d0690884d3219a13a4e'
def load_penalty_after_authorization():
 assert sha(PENALTY_PATH)==PENALTY_SHA;spec=importlib.util.spec_from_file_location('verified153_penalty155',PENALTY_PATH);api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api);return api
def new_parameters_after_authorization(source,log_upper,api):
 c=source['old_core'];e=F(c['epsilon']);old=F(c['penalty_coefficient']);U=F(source['old_H_interval'][1]);assert U>0 and F(c['target110'])==F(11,100000);mass=source['mass_qualification']
 if e==F(1,64):assert source['new_rule']=='153_EPS64_uniform_no_H32_mass' and old==256 and mass is None;L=F(0)
 elif c['curvature_rule_class']=='075_full_frame':assert e==F(1,128) and old==1024 and mass is None;L=F(0)
 else:
  assert e==F(1,128) and c['curvature_rule_class']=='129_anchor_fixed_Jensen' and mass is not None;assert mass['lower_epsilon']=='1/64' and mass['upper_epsilon']=='1/128' and mass['point']==c['point'] and F(mass['upper'])==U;L=F(mass['lower']);assert L==F(21,200000) and L<U and old==1024-768*L/U and mass['point_in_support_closedbox'] and mass['no_new_raw_lower_lookup'] and mass['global_H128_floor_not_used_as_H64']
 caps=api.curvature_caps(e,L,U,old);q=caps['new_q_upper'];assert 0<q<=old and q<=caps['uniform_time_cap']
 if e==F(1,64):assert q==F(512,3)
 if e==F(1,128) and mass is None:assert q==F(2048,3)
 # No legacy129 lower bound (256 or A/4 or A/6) is imposed on newq.
 return {'canonical_id':source['canonical_id'],'epsilon':c['epsilon'],'point':c['point'],'target110':'11/100000','target115':'23/200000','Delta110_lower':c['Delta110_lower'],'Delta115_lower':str(F(c['Delta110_lower'])-log_upper),'loggradient_directed':c['loggradient_directed'],'q_old':str(old),'q_new':str(q),'q_new_rule':'153_time_shape_uniform_or_original_g105_mass','caps153':{k:(str(v) if v is not None else None) for k,v in caps.items()},'mass_qualification':mass,'raw_record_sha256':source['raw_record_sha256'],'source_penalty_sha256':PENALTY_SHA}
def one_new16corner_calculation(parameters,leaf,log_upper):
 assert parameters['canonical_id']==leaf['canonical_id'] and parameters['epsilon']==leaf['epsilon'];point=list(map(F,parameters['point']));delta=F(parameters['Delta110_lower']);gamma=[tuple(map(F,g)) for g in parameters['loggradient_directed']];q=F(parameters['q_new']);lo=list(map(F,leaf['box_lower']));hi=list(map(F,leaf['box_upper']));assert q>0 and len(gamma)==4 and all(a<=b for a,b in gamma)
 rows=[]
 for bits in product((0,1),repeat=4):
  v=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)];rows.append((delta+sum(min(a*x,b*x) for x,(a,b) in zip(v,gamma))-q*sum(x*x for x in v),bits))
 margin,bits=min(rows);old=F(leaf['saved_g110_margin']);assert margin>=old>0;shifted=margin-log_upper;assert F(parameters['Delta115_lower'])==delta-log_upper
 return {'leaf_index':leaf['leaf_index'],'node':leaf['node'],'canonical_id':leaf['canonical_id'],'epsilon':leaf['epsilon'],'q_new_rule':parameters['q_new_rule'],'q_old':parameters['q_old'],'q_new':parameters['q_new'],'margin110_new':str(margin),'margin115':str(shifted),'saved_old_margin110':leaf['saved_g110_margin'],'margin110_new_ge_old_positive':True,'worst_corner_bits':list(bits),'checked_points':16,'passed115':shifted>0,'full110_leaf_preserved':True,'one_anchor_all16corners':True}
