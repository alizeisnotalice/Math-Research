#!/usr/bin/env python3
"""Independent saved146 integer audit. No Joint/producer/Arb imports or root calculations."""
import sys
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from fractions import Fraction as F
from math import comb
from concurrent.futures import ProcessPoolExecutor,as_completed
import gzip,hashlib,json,multiprocessing,resource,time,traceback
R=Path(__file__).resolve().parent
checks=0

def ck(v,msg):
 global checks
 checks+=1
 if not v:raise ArithmeticError(msg)
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(v):return json.dumps(v,sort_keys=True,separators=(',',':')).encode()
def save(n,d):
 p=R/n;t=p.with_name(p.name+'.own-tmp')
 with t.open('w')as f:json.dump(d,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def fl(q,S):return q.numerator*S//q.denominator
def ce(q,S):return-((-q.numerator*S)//q.denominator)
def unrank(index):
 # Independent lexicographic combinatorial unranking; no islice or producer iterator.
 ck(0<=index<comb(24,8),'valid prefix rank');out=[];minimum=0
 for remaining in range(8,0,-1):
  for v in range(minimum,17):
   block=comb(16-v+remaining-1,remaining-1)
   if index<block:out.append(v);minimum=v;break
   index-=block
  else:raise ArithmeticError('unranking exhausted')
 ck(index==0,'unrank exact remainder');return out

def advance(state):
 for k in range(7,-1,-1):
  if state[k]<16:
   z=state[k]+1
   return state[:k]+[z]*(8-k)
 return None

def probability(state,powers):
 cats=[state[0]]+[state[j]-state[j-1]for j in range(1,8)]+[32-state[-1]]
 ck(min(cats)>=0 and sum(cats)==32,'true category vector')
 remaining=32;coefficient=1;out=1
 for j,c in enumerate(cats):
  coefficient*=comb(remaining,c);remaining-=c;out*=powers[j][c]
 ck(remaining==0,'sequential binomial coefficient complete');return coefficient*out

def worker(payload):
 global checks
 checks=0;job,table,rootpath=payload['header'],payload['table'],Path(payload['root_path']);start=time.monotonic();c0=time.process_time()
 ck(os.getpriority(os.PRIO_PROCESS,0)>=10,'worker nice10');ck(sha(rootpath)==job['root_file_sha256'],'unchanged saved new root chunk')
 S=1<<192;S6=S**6;common=list(map(int,table['common_W_units']));groups=[[[int(z)for z in pair]for pair in row]for row in table['group_W_units']];powers=[[[int(q)**n for n in range(33)]for q in row]for row in table['category_units']]
 with gzip.open(rootpath,'rt')as f:records=json.load(f)
 ck(len(records)==job['states']==job['end_exclusive']-job['start'],'exact chunk count');state=unrank(job['start']);first=state.copy();last=None;digest=hashlib.sha256();countK=[0]*17;mass=[0,0];low=upper=0
 for record in records:
  ck(state is not None and len(record)==2,'complete state/root pair')
  digest.update(bytes(state));wl=common[0]+sum(groups[j][n][0]for j,n in enumerate(state));wu=common[1]+sum(groups[j][n][1]for j,n in enumerate(state));vl,vu=map(int,record)
  ck(0<wl<=wu and 0<vl<=vu,'positive W and root integer interval');ck(vl**5*wu<=S6<=vu**5*wl,'saved exact fifth-root cross products')
  pl,pu=probability(state,powers[0]),probability(state,powers[1]);low+=pu*vl;upper+=pl*vu;mass[0]+=pl;mass[1]+=pu;countK[state[-1]]+=1;last=state;state=advance(state)
 ck(first==job['first']and last==job['last'],'actual chunk endpoint identities');ck(digest.hexdigest()==job['state_bytes_sha256'],'independent actual chunk state hash')
 ck(countK==job['prefix_count_by_K'],'actual per-K count');ck(mass==list(map(int,job['probability_numerators_by_law'])),'both full-law chunk masses')
 ck(low==int(job['lower_at_pU_numerator'])and upper==int(job['upper_at_pL_numerator']),'independent positive weighted chunk sums')
 if job['end_exclusive']<735471:ck(state==unrank(job['end_exclusive']),'precise chunk next index')
 else:ck(state is None,'entire prefix exhausted')
 ck(sha(rootpath)==job['root_file_sha256'],'root source unchanged after read')
 return{'job':job['job'],'states':len(records),'first':first,'last':last,'state_bytes_sha256':digest.hexdigest(),'countK':countK,'mass':list(map(str,mass)),'low':str(low),'upper':str(upper),'root_checks':len(records),'assertions':checks,'CPU_seconds':time.process_time()-c0,'wall_seconds':time.monotonic()-start}

def main():
 start=time.monotonic();c0=time.process_time();spec=rd(R/'inputs.json');reg=rd(R/'registration.json');sources=rd(R/'source_manifest.json')['sources'];status='FAILED_STOP_NO_RETRY';failure=None;parts=[]
 ck(reg['checker_sha256']==sha(__file__)and reg['inputs_sha256']==sha(R/'inputs.json')and reg['source_manifest_sha256']==sha(R/'source_manifest.json'),'registered exact new audit')
 for row in sources:ck(sha(row['path'])==row['sha256'],'pre source '+row['path'])
 ck(not(R/'started.json').exists(),'sole independent audit, no automatic restart');n=os.getpriority(os.PRIO_PROCESS,0)
 if n<10:os.nice(10-n)
 with(R/'started.json').open('x')as f:json.dump({'pid':os.getpid(),'unix_time':time.time(),'registration_sha256':sha(R/'registration.json'),'checker_sha256':sha(__file__)},f);f.write('\n');f.flush();os.fsync(f.fileno())
 try:
  I=rd(spec['author_inputs_path']);cert=rd(spec['certificate_path']);tables=rd(spec['tables_path']);jobs=rd(spec['jobs_path']);term=rd(spec['main_terminal_path']);ex=rd(spec['author_exit_path']);old=rd(spec['reference128_input_path']);r144=rd(spec['reference144_result_path']);o144=rd(spec['reference144_input_path']);r130=rd(spec['reference130_result_path']);oldcert=rd(spec['reference128_certificate_path'])
  ck(term['exit_code']==0 and term['restarts']==0 and term['source_sha256']==spec['author_driver_sha256']and term['release_sha256']==sha(spec['author_release_path']),'actual primary original exit/source/release')
  ck(ex['certificate_sha256']==sha(spec['certificate_path'])and ex['workers_completed']==15,'original terminal certificate links')
  ck(cert['status']=='COMPLETE_ONE_NEW_RAW_H_AND_TRUE_COST_CERTIFICATE'and cert['candidate']==I['candidate']==r144['candidate'],'actual new candidate')
  ck(I['raw_model_input']==o144,'raw optimistic numerical model source')
  cand=I['candidate'];den=int(cand['denominator']);nums=cand['numerators'];ck(den==2**44 and len(nums)==49 and min(nums)>=0 and nums[0]>0 and sum(nums)==den,'exact44-bit nonnegative simplex')
  x=[F(n,den)for n in nums];raw=I['raw_model_input'];tb=[F(raw['beta0'])]+[F(v)for row in raw['beta_by_h']for v in row];caps=[F(old['beta0'])]+[F(v)for row in old['beta_by_h']for v in row];a=[z/b for z,b in zip(x,tb)]
  ck(len(a)==len(caps)==len(tb)==49 and list(map(str,a))==I['raw_coefficients_exact'],'all exact unquantized raw coefficients')
  ck(list(map(F,I['valid_cap_vector_exact']))==caps,'all49 legal source caps');ck(I['changed_index']==32 and F(I['x_changed_exact'])==x[32]>0,'unique changed feature')
  bL=F(I['bL_definition_only']);b=F(I['current_valid_beta2_upper']);ck(bL==F(r130['moment_interval'][0])==tb[32]and b==min(caps[32],F(r130['moment_interval'][1]))==caps[32],'legal b cap and lower-only bL')
  ck(all(tb[j]==caps[j]for j in range(49)if j!=32),'other legal caps unchanged');C=sum(aa*bb for aa,bb in zip(a,caps));ck(C==F(I['current_cost_exact'])==1-x[32]+x[32]*b/bL==F(cert['current_true_cost_upper_exact'])>0,'independent exact49-term legal cost')
  ck(cert['raw_a_quantized']is False and raw['Phi_table']==old['Phi_table']and raw['p_intervals']==old['p_intervals']and raw['h_values']==old['h_values'],'same Phi/law/h, no a rounding')
  S=1<<192;D=S**33;S6=S**6;phi=old['Phi_table'];common=[fl(a[0],S),ce(a[0],S)];ck(common[0]>0,'positive intercept lower')
  aa=[a[1+6*j:1+6*(j+1)]for j in range(8)];group=[]
  for j in range(8):
   rows=[]
   for n in range(33):
    value=sum(aa[j][k]*phi[n][k+1]for k in range(6));rows.append([fl(value,S),ce(value,S)])
   group.append(rows)
  ps=[];qs=[]
  for side in range(2):
   p=[]
   for pair in raw['p_intervals']:
    value=F(pair[side])*S;ck(value.denominator==1,'exact192-grid CDF endpoint');p.append(int(value))
   q=[p[0]]+[p[j]-p[j-1]for j in range(1,8)]+[S-p[-1]];ck(min(q)>0 and sum(q)==S,'whole honest endpoint law');ps.append(p);qs.append(q)
  ck(all(l<=u for l,u in zip(ps[0],ps[1])),'componentwise CDF coupling')
  expected={'scale_power2':192,'common_W_units':list(map(str,common)),'group_W_units':[[list(map(str,pair))for pair in row]for row in group],'p_units':[[str(v)for v in row]for row in ps],'category_units':[[str(v)for v in row]for row in qs],'probability_denominator_power2':6144,'weighted_sum_denominator_power2':6336,'a0':str(a[0]),'a_by_h':[[str(v)for v in row]for row in aa]}
  ck(tables==expected,'all264 rational W groups/intercept/a/CDF tables independently rebuilt');tablehash=hashlib.sha256(canon(expected)).hexdigest()
  ck(len(jobs)==len(cert['chunks'])==15 and jobs[0]['start']==0 and jobs[-1]['end_exclusive']==735471,'exact15 chunks full prefix')
  for j,(job,h)in enumerate(zip(jobs,cert['chunks'])):
   ck(job['job']==h['job']==j and 0<=job['start']<job['end_exclusive']<=735471 and job['states']==job['end_exclusive']-job['start'],'proper independent chunk range')
   ck(all(job[k]==h[k]for k in job),'actual header/planned job identity');ck(h['integer_tables_sha256']==tablehash,'actual shared raw W table hash')
   if j:ck(jobs[j-1]['end_exclusive']==job['start'],'no gap or overlapping new prefix')
  payloads=[{'header':h,'table':expected,'root_path':str(Path(spec['author_directory'])/h['root_file'])}for h in cert['chunks']]
  with ProcessPoolExecutor(max_workers=15,mp_context=multiprocessing.get_context('spawn'))as pool:
   futures=[pool.submit(worker,p)for p in payloads]
   for f in as_completed(futures):
    parts.append(f.result());save('progress.json',{'phase':'independent_saved_integer_chunks','completed_chunks':len(parts),'prefix_root_checks':sum(p['root_checks']for p in parts)})
  parts.sort(key=lambda p:p['job']);ck(sum(p['root_checks']for p in parts)==735471,'every prefix root inequality checked')
  # Independent enumeration-only global identity, no further roots or objective calls.
  digest=hashlib.sha256();state=[0]*8;countK=[0]*17;chunkdig=hashlib.sha256();idx=0;ji=0
  while state is not None:
   bstate=bytes(state);digest.update(bstate);chunkdig.update(bstate);countK[state[-1]]+=1;idx+=1
   if idx==jobs[ji]['end_exclusive']:
    ck(chunkdig.hexdigest()==parts[ji]['state_bytes_sha256'],'independent contiguous whole-state stream hash');ji+=1;chunkdig=hashlib.sha256()
   state=advance(state)
  ck(idx==735471 and ji==15 and countK==[comb(k+7,7)for k in range(17)]==cert['prefix_count_by_K'],'complete independent prefix/per-K identity')
  ck(digest.hexdigest()==I['expected_state_hash']==cert['state_bytes_sha256'],'complete new state hash')
  ck(all(sum(p['countK'][k]for p in parts)==countK[k]for k in range(17)),'sum independent worker per-K counts')
  mass=[sum(int(p['mass'][s])for p in parts)for s in range(2)];totalL=sum(int(p['low'])for p in parts);totalU=sum(int(p['upper'])for p in parts)
  ck(mass==list(map(int,cert['prefix_probability_numerators_by_law'])),'complete two prefix-law masses')
  tailMass=[0,0];tailU=0;ck(len(cert['tail_records'])==16,'all16 true tails')
  for K,row in zip(range(17,33),cert['tail_records']):
   ck(row['K']==K,'tail K order');wl=common[0]+group[-1][K][0];wu=common[1]+group[-1][K][1];ck([wl,wu]==list(map(int,row['W_units'])),'actual raw tail W enclosure')
   vl,vu=map(int,row['v_units']);ck(0<vl<=vu and vl**5*wu<=S6<=vu**5*wl,'exact saved tail root inequality')
   pp=[comb(32,K)*ps[s][-1]**K*(S-ps[s][-1])**(32-K)for s in range(2)];ck(pp==list(map(int,row['probability_numerators_by_law'])),'exact two binomial tail laws');tailMass=[tailMass[s]+pp[s]for s in range(2)];tailU+=pp[0]*vu
  ck(tailMass==list(map(int,cert['tail_probability_numerators_by_law']))and all(mass[s]+tailMass[s]==S**32 for s in range(2)),'both entire endpoint laws mass1')
  HL,HU=F(totalL,D),F(totalU+tailU,D);ck(0<HL<=HU and[HL,HU]==list(map(F,cert['H_interval'])),'complete independent positive H interval')
  L5,U5=HL**5,HU**5;ck([L5,U5]==list(map(F,cert['H_power5_interval'])),'exact H fifth powers')
  Mnew=C*U5;ck(Mnew==F(cert['current_valid_M6_power6_upper']),'legal cost restores actual valid M upper budget')
  theta=F(8,9)**6*24576**5/65000**6;old6=F(oldcert['M_power6_upper']);ck(theta==F(I['threshold_g0_power6_exact'])==F(old['threshold_M_power6_exact']),'full accurate g0 threshold')
  ck(cert['current_cap_beats_g0_threshold']==(Mnew<theta)and cert['current_cap_improves_saved128']==(Mnew<old6),'accurate two cap-budget comparisons')
  def threshold(t,q):return bL*(t/q-(1-x[32]))/x[32]
  thresholds={}
  for label,t,field in [('g0',theta,'true_beta2_upper_threshold_for_strict_g0'),('saved128_improvement',old6,'true_beta2_upper_threshold_for_strict_saved128_improvement')]:
   suff,necessary=threshold(t,U5),threshold(t,L5);ck(suff<=necessary,'correct sufficient/necessary threshold order');ck(F(cert[field]['sufficient_exact'])==suff and F(cert[field]['necessary_for_this_raw_candidate_exact'])==necessary,'independent exact legal beta2 thresholds')
   ck((C*U5<t)==(b<suff)and(C*L5<t)==(b<necessary),'strict beta2 rearrangement directions');thresholds[label]={'sufficient_exact':str(suff),'necessary_for_fixed_raw_cap_budget_exact':str(necessary),'current_legal_b_below_sufficient':b<suff,'current_legal_b_below_necessary':b<necessary}
  ck(F(cert['full_I_upper_power5_from123'])==F(3,16)**6*Mnew/384**5,'complete I budget cost once')
  ck(cert['logical_root_requests']==735487 and cert['prefix_states']==735471 and cert['probability_denominator_power2']==6144 and cert['weighted_sum_denominator_power2']==6336,'actual fixed certificate accounting')
  ck(cert['both_complete_endpoint_laws_unit_mass']is True and cert['tail_discarded']is False and cert['probability_renormalized']is False,'no omitted or renormalized probability')
  ck(cert['optimizer_kernel_erfc_native_calls']==0 and cert['old_baseline_replays']==0 and cert['new_C6_bound']is False and cert['global_g0_proved']is False,'no false native/baseline/global claims')
  for row in sources:ck(sha(row['path'])==row['sha256'],'post source '+row['path'])
  status='PASS_COMPLETE_SAVED_NEW146_RAW_H_AND_LEGAL_COST'
  result={'status':status,'assertions':checks+sum(p['assertions']for p in parts),'prefix_root_inequalities_checked':735471,'tail_root_inequalities_checked':16,'total_saved_root_inequalities':735487,'chunks':15,'independent_workers_max':15,'source_files_pre_post':len(sources),'independent_probability_method':'sequential-binomial multinomial coefficients; combinatorial unranking and next-state advancement, no producer imports','H_interval':[str(HL),str(HU)],'H_power5_interval':[str(L5),str(U5)],'C_legal_cost_exact':str(C),'current_legal_M_power6_upper':str(Mnew),'old128_M_power6_upper':str(old6),'g0_threshold_exact':str(theta),'current_cap_beats_g0':Mnew<theta,'current_cap_improves_saved128':Mnew<old6,'beta2_thresholds':thresholds,'full_I_upper_power5':str(F(3,16)**6*Mnew/384**5),'both_complete_laws_mass1':True,'state_hash':digest.hexdigest(),'actual_M_or_C6_lower_bound':False,'global_g_or_C6_proved':False,'new_Arb_roots_or_old_baseline_replays':0,'worker_CPU_sum':sum(p['CPU_seconds']for p in parts),'worker_summaries':parts,'wall_before_receipt':time.monotonic()-start,'collector_CPU_before_receipt':time.process_time()-c0}
  save('review_results.json',result)
 except BaseException as e:
  failure=repr(e);save('failure.json',{'failure':failure,'traceback':traceback.format_exc(),'completed_saved_chunks':len(parts),'automatic_retry':False});raise
 finally:save('run_receipt.json',{'status':status,'pid':os.getpid(),'failure':failure,'wall_before_receipt':time.monotonic()-start,'collector_CPU_seconds':time.process_time()-c0,'completed_saved_chunks':len(parts),'worker_CPU_sum':sum(p['CPU_seconds']for p in parts),'actual_max_science_workers':15,'priority':os.getpriority(os.PRIO_PROCESS,0),'time_or_CPU_ceiling':None,'automatic_retries':0,'new_Arb_or_native_calls':0})
 print(json.dumps({'status':status,'assertions':checks+sum(p['assertions']for p in parts),'saved_root_inequalities':735487,'current_cap_beats_g0':Mnew<theta,'current_cap_improves_saved128':Mnew<old6}),flush=True)
if __name__=='__main__':main()
