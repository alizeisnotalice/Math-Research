"""C161 only NEW scrambled Sobol/Exp6 estimates; import has no scientific side effects."""
def evaluate_job(job,cone,packet,checkpoint,counters,sample_state):
 import hashlib,math,resource,time
 from fractions import Fraction as F
 import numpy as np
 import scipy
 from scipy.special import erfcx
 from scipy.stats import qmc
 from common import need
 need(np.__version__==packet['versions']['numpy']and scipy.__version__==packet['versions']['scipy'],'PREPARED_NUMPY_SCIPY_VERSIONS')
 start=time.monotonic();cpu=time.process_time();counters.update({'sobol_jobs_started':0,'sample_points':0,'uniform_coordinates':0,'log1p_scalar_evaluations':0,'kernel_scalar_evaluations':0,'erfcx_scalar_evaluations':0,'exp_scalar_evaluations':0,'positive_prefix_updates':0,'positive_pair_prefix_updates':0})
 h=float(F(packet['h']));s=math.sqrt(h);B=np.array([[float(F(v))for v in row]for row in cone['nodes_by_body']],dtype=np.float64)
 alpha=float(F(cone['alpha']));pL,pU=map(F,packet['p_source_rational_enclosure']);p=float((pL+pU)/2)
 need(B.shape==(32,6)and np.all(np.isfinite(B))and np.all(B>=0),'NUMERICAL_SOURCE_SHAPE')
 need(packet['sobol_m']==14 and packet['chunk_points']==16384 and packet['sobol_bits']==53,'FIXED_ONE_CHUNK')
 checkpoint('rng_intent.json',{'job':job,'method':'scipy.stats.qmc.Sobol(d=6,scramble=True,bits=53,seed=seed).random_base2(m=14)','versions':packet['versions']})
 counters['sobol_jobs_started']+=1
 sample_state.update({'phase':'SCRAMBLE_CONSTRUCTOR_INFLIGHT','sobol_constructor_attempted':True})
 engine=qmc.Sobol(d=6,scramble=True,bits=53,seed=job['seed'])
 sample_state.update({'phase':'SOBOL_REQUEST_INFLIGHT','sobol_request_attempted':True,'sobol_request_completed':False,'actual_sample_points':None,'actual_sample_points_known':False})
 counters['sample_points']=None;counters['uniform_coordinates']=None
 U=engine.random_base2(m=14)
 sample_state.update({'phase':'SOBOL_RETURNED_UNVALIDATED','sobol_request_completed':True,'actual_sample_points':None,'actual_sample_points_known':False})
 need(U.shape==(16384,6)and np.all(np.isfinite(U))and np.all(U>=0)and np.all(U<1),'TRUE_UNIT_SOBOL_SAMPLES')
 counters['sample_points']=len(U);counters['uniform_coordinates']=int(U.size)
 sample_state.update({'phase':'COMPLETE_VALIDATED_16384_POINT_SAMPLE','actual_sample_points':len(U),'actual_sample_points_known':True})
 uhash=hashlib.sha256(np.asarray(U,dtype='<f8',order='C').tobytes(order='C')).hexdigest()
 checkpoint('samples_receipt.json',{'job_id':job['job_id'],'seed':job['seed'],'U_sha256':uhash,'U_hash_encoding':'little-endian IEEE754 float64 C row-major','shape':list(U.shape),'min_U':float(U.min()),'max_U':float(U.max()),'zero_coordinates':int(np.count_nonzero(U==0)),'samples_are_new':True})
 with np.errstate(divide='raise',invalid='raise',over='raise',under='ignore'):
  counters['log1p_scalar_evaluations']+=int(U.size);E=-np.log1p(-U)
  need(np.all(np.isfinite(E))and np.all(E>=0),'FINITE_EXP6_TRANSFORM')
  prefix=np.zeros(len(U),dtype=np.float64);pairs=np.zeros(len(U),dtype=np.float64);kmin=math.inf;kmax=-math.inf;zero_k=0
  # Same E row drives ALL32 forms. No pairwise/body/cone pseudo independence.
  for b in B:
   a=E@b;u=a/(2*s);left=u-s;right=u+s
   counters['erfcx_scalar_evaluations']+=2*len(U)
   difference=erfcx(left)-erfcx(right)
   counters['exp_scalar_evaluations']+=len(U);counters['kernel_scalar_evaluations']+=len(U)
   k=0.5*np.exp(-(left*left))*difference
   need(np.all(np.isfinite(k))and np.all(k>=0),'NUMERICAL_KERNEL_NONNEGATIVE_SANITY_NO_CLIPPING')
   need(np.all(k<=p*(1+1e-12)+1e-15),'NUMERICAL_KERNEL_RANGE_SANITY_NOT_CERTIFICATION')
   kmin=min(kmin,float(k.min()));kmax=max(kmax,float(k.max()));zero_k+=int(np.count_nonzero(k==0))
   # e2 of32 values:32 multiplications/additions, all nonnegative; no subtractive square.
   counters['positive_pair_prefix_updates']+=len(U);pairs+=prefix*k
   counters['positive_prefix_updates']+=len(U);prefix+=k
  need(np.all(np.isfinite(pairs))and np.all(pairs>=0),'FINITE_NONNEGATIVE_PAIR_SAMPLE')
  mean=float(np.mean(pairs,dtype=np.float64));mean_square=float(np.mean(pairs*pairs,dtype=np.float64))
 return {'status':'COMPLETE_RQMC_CONE_DIAGNOSTIC','job':job,'U_sha256':uhash,'alpha_exact':cone['alpha'],'alpha_float_hex':alpha.hex(),'cone_mean_e2_float':mean,'weighted_cone_mean_float':alpha*mean,'sample_second_raw_moment_e2_float':mean_square,'sample_min_e2_float':float(pairs.min()),'sample_max_e2_float':float(pairs.max()),'sample_kernel_min_float':kmin,'sample_kernel_max_float':kmax,'zero_kernel_values':zero_k,'E_min_float':float(E.min()),'E_max_float':float(E.max()),'h_float_hex':h.hex(),'numerical_p_reference_float_hex':p.hex(),'actual_counts':dict(counters),'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu,'worker_maxrss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'worker_maxrss_units':'bytes on macOS; KiB on Linux','strict_mu2_lower':False,'strict_mu2_upper':False,'new_Cn_bound':False,'density_generations':0,'old_grid_MC_replays':0}
