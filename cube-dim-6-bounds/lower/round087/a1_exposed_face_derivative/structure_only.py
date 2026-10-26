#!/usr/bin/env python3
"""Two structural germs only after DP cap: zero section/union nodes."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,importlib.util,time,signal,resource,sys,os
sys.dont_write_bytecode=True
H=Path(__file__).resolve().parent;OUT=H/'structure_results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21));signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
S=H/'results.json';IF=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/lower/round008/spacing_chamber/chamber_interface.json');AR=IF.parents[2]/'round001_rectangular_certificate.py'
# IF.parents[2] is lower/; explicit known source path avoids private inference.
AR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/lower/round001_rectangular_certificate.py')
s=json.loads(S.read_text());a=json.loads(IF.read_text());assert s['status']=='stopped_exception' and s['new_DP_nodes']==1001
spec=importlib.util.spec_from_file_location('accepted_field',AR);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
A.D=10**80;A.RU,A.RV=A.floor_root(2*A.D**6,6),A.floor_root(3*A.D**6,6);A.UL,A.UH=Q(A.RU,A.D),Q(A.RU+1,A.D);A.VL,A.VH=Q(A.RV,A.D),Q(A.RV+1,A.D);A.bounds.cache_clear()
F,ZERO=A.F,A.F(0)
def field(z):
 x=[Q(0)]*36
 for ij,c in z.items():i,j=map(int,ij.split(','));x[6*i+j]=Q(c)
 return F(a=x)
def enc(x):return {f'{i},{j}':str(x.a[6*i+j]) for i in range(6) for j in range(6) if x.a[6*i+j]}
def view(x):
 lo,hi=A.bounds(x.a);return {'coefficients':enc(x),'exact_interval':[str(lo),str(hi)],'display':float((lo+hi)/2)}
funcs=[{'c':field(z['c']),'m':Q(z['m'])} for z in s['normal_functions']];t=field(s['a1_current']['coefficients']);moving={z['raw_index']:i for i,f in enumerate(s['normal_functions']) for z in f['labels'] if 'raw_index' in z}
H5=(field(s['fixed_tail_spacing']['coefficients'])/2)**5
out={'status':'finished_structure_only','pid':os.getpid(),'new_configurations':2,'new_DP_nodes':0,'new_5D_union_calls':0,'static_6D_union_calls':0,'prior_DP_admitted_nodes':1000,'prior_DP_counter_including_rejected_request':1001,'new_C6_bound':False,'sides':{}}
for side,label in [(-1,'left'),(1,'right')]:
 keys=[(f['c']+f['m']*t,side*f['m']) for f in funcs];sel={};used=set()
 for k,ell in enumerate(a['box_lengths']):
  lo=max((0,) if ell[0]%2 else (0,moving[k]),key=lambda i:keys[i]);hi=min((1,moving[k]) if ell[0]%2 else (1,),key=lambda i:keys[i])
  if keys[lo]<keys[hi]:sel[k]=(lo,hi);used.update((lo,hi))
 order=sorted(used,key=lambda i:keys[i]);allorder=sorted(range(len(funcs)),key=lambda i:keys[i]);pos={x:i for i,x in enumerate(order)};ties=[];block=[allorder[0]]
 for i,j in zip(allorder,allorder[1:]):
  if keys[i][0]==keys[j][0]:block.append(j)
  else:
   if len(block)>1:ties.append(block)
   block=[j]
 if len(block)>1:ties.append(block)
 slabs=[];neg=Q(0);positive=Q(0)
 for lo,hi in zip(order,order[1:]):
  c=funcs[hi]['m']-funcs[lo]['m']
  if not c:continue
  active=[k for k,(p,q) in sel.items() if pos[p]<=pos[lo] and pos[hi]<=pos[q]]
  slabs.append({'ids':[lo,hi],'width_slope':str(c),'active_raw_indices':active,'zero_baseline_width':keys[lo][0]==keys[hi][0],'5D_union_area':'NOT_COMPUTED'})
  neg+=min(c,0);positive+=max(c,0)
 out['sides'][label]={'positive_boxes':len(sel),'selections':{str(k):list(z) for k,z in sel.items()},'effective_order':order,'all_function_order':allorder,'all_tie_groups':ties,'required_weighted_slabs':slabs,'trivial_derivative_lower':view(64*neg*H5),'trivial_derivative_upper':view(64*positive*H5),'trivial_bound_basis':'0<=each complete5D section union area<=half-tail-cell-volume; no overlap assumption or new section integration'}
out['minmax_selections_unchanged_across_current']=out['sides']['left']['selections']==out['sides']['right']['selections']
out.update({'CPU_seconds':time.process_time()-CPU,'combined_CPU_seconds':s['CPU_seconds']+time.process_time()-CPU,'wall_seconds_before_save':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S,IF,AR,Path(__file__)]}})
with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps({k:out[k] for k in ['status','pid','new_configurations','new_DP_nodes','CPU_seconds','combined_CPU_seconds','wall_seconds_before_save']}));signal.alarm(0)
