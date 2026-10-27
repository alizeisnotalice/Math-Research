#!/usr/bin/env python3
"""Terra prepare-only: exact input identity/finite registration, NO scientific imports."""
import sys
sys.dont_write_bytecode=True
from fractions import Fraction as F
from math import prod
import hashlib,importlib.metadata,itertools,json,platform
from common import P,need,read,sha,write_new,check_sources,source_freeze

def prepare():
 need(sys.argv[1:]==['--prepare-only'],'PREPARE_ONLY_FLAG_REQUIRED')
 release=read(P/'prepare_release.json')
 need(release.get('kind')=='C161_SOURCE_PREPARE_ONLY'and release.get('main_approved')is True and release.get('source_freeze')==source_freeze(),'SEPARATE_PREPARE_ONLY_GATE')
 need(not(P/'prepare_started.json').exists(),'NO_PREPARE_RESTART_OR_OVERWRITE')
 write_new(P/'prepare_started.json',{'status':'PREPARE_ONLY_NO_RNG_KERNEL_OR_OLD_BASELINE','new_samples':0})
 pins=read(P/'source_pins.json');check_sources(pins)
 spec=read(pins['sources']['130_inputs']['path']);table=read(pins['sources']['130_cone_table']['path']);old=read(pins['sources']['130_results']['path'])
 need(F(spec['h'])==F(1,128)and spec['cone_count']==21 and spec['body_count']==32 and spec['pair_count']==496,'PINNED_TRUE_130_MODEL')
 cones=table['cones'];need(len(cones)==21 and [c['id']for c in cones]==list(range(21)),'ALL21_SOURCE_CONES')
 V=[(1,)+tail for tail in itertools.product((-1,1),repeat=5)]
 for c in cones:
  need(len(c['rays'])==len(c['rates'])==6 and all(len(r)==6 for r in c['rays']),'EXACT_RAW_CONE_SHAPE')
  need(all(type(r)is int and r>0 for r in c['rates']),'POSITIVE_SOURCE_RATES')
  need(F(c['alpha'])==F(c['determinant'],prod(c['rates']))>0,'EXACT_ALPHA_FROM_SOURCE')
  need(len(c['nodes_by_body'])==32 and all(len(row)==6 for row in c['nodes_by_body']),'ALL32_BODY_ROWS')
  for row,v in zip(c['nodes_by_body'],V):
   values=list(map(F,row));expected=[F(abs(sum(v[k]*ray[k]for k in range(6))),rate)for ray,rate in zip(c['rays'],c['rates'])]
   need(values==expected and any(values)and all(0<=x<F(1,2)for x in values),'ORIGINAL32_SPECTRAL_NODES')
 need(46080*sum((F(c['alpha'])for c in cones),F(0))==F(spec['R0']),'FULL_CONE_WEIGHT_NORMALIZATION')
 kappa=list(map(F,old['pi_factor_interval']));need(0<kappa[0]<=kappa[1],'PINNED_SYMBOLIC_46080_OVER_PI3_ENCLOSURE')
 p=list(map(F,spec['p_interval']));need(0<p[0]<=p[1]<1,'PINNED_P_RANGE')
 source_identity=sha(pins['sources']['130_cone_table']['path']);jobs=[];seeds=set()
 for rep in range(8):
  for cone in range(21):
   label='C161_RQMC_v1|rep='+str(rep)+'|cone='+str(cone)+'|table='+source_identity
   digest=hashlib.sha256(label.encode()).hexdigest();seed=int(digest[:8],16)
   need(seed not in seeds,'DISTINCT_SCRAMBLE_SEEDS');seeds.add(seed)
   jobs.append({'job_id':'r'+str(rep)+'_c'+str(cone).zfill(2),'replicate':rep,'cone_id':cone,'seed':seed,'seed_label':label,'seed_label_sha256':digest})
 versions={name:importlib.metadata.version(name)for name in ('numpy','scipy')}
 packet={'status':'PREPARED_FINITE_RQMC_DIAGNOSTIC_UNRUN','h':spec['h'],'kappa_symbolic':'46080/pi^3','kappa_source_rational_enclosure':old['pi_factor_interval'],'kappa_float_rule':'float(exact rational enclosure midpoint), diagnostic only','p_source_rational_enclosure':spec['p_interval'],'cones':cones,'jobs':jobs,'replicates':8,'cone_count':21,'body_count':32,'sobol_dimension':6,'sobol_m':14,'sobol_bits':53,'chunk_points':16384,'versions':versions,'python_version':platform.python_version(),'source_pins_sha256':sha(P/'source_pins.json'),'source_table_sha256':source_identity,'source_freeze':source_freeze(),'new_density_dependency':False,'old_kernel_or_grid_replays':0}
 write_new(P/'case_packet.json',packet)
 plan=read(P/'plan.json')
 registration={'status':'PREPARED_UNRUN_PENDING_SCIENCE_DUAL_GATE','plan_sha256':sha(P/'plan.json'),'case_packet_sha256':sha(P/'case_packet.json'),'source_pins_sha256':sha(P/'source_pins.json'),'source_freeze':packet['source_freeze'],'jobs':168,'replicates':8,'finite_counts':plan['finite_counts'],'max_workers':15,'nice_minimum':19,'time_limit':None,'automatic_retries':0,'versions':versions,'python_version':packet['python_version'],'numerical_estimate_only':True,'prepare_rng_or_kernel_calls':0}
 write_new(P/'registration.json',registration)
 write_new(P/'prepare_receipt.json',{'status':'COMPLETE_PREPARE_ONLY','registration_sha256':sha(P/'registration.json'),'case_packet_sha256':sha(P/'case_packet.json'),'jobs':168,'new_samples':0,'rng_kernel_special_function_calls':0})
 print(json.dumps({'status':'COMPLETE_PREPARE_ONLY','jobs':168,'registration_sha256':sha(P/'registration.json')}))
if __name__=='__main__':prepare()
