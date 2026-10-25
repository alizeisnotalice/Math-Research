#!/usr/bin/env python3
"""Exact full-group Green plus at most two H passes; stop at first true refutation."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from time import perf_counter,time
import hashlib,importlib.util,itertools,json,math,os,signal
import flint
from flint import arb,fmpq,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
SOURCE=PROJECT/'upper/round007/g_time_moments/run_moments.py';FAN=PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
index=int(sys.argv[1])if len(sys.argv)>1 else 1
assert index in (1,2)
PFILE=HERE/f'point{index}.json';point=json.loads(PFILE.read_text());gr=json.loads((HERE/'gradient_results.json').read_text())
remaining=1200-(time()-gr['batch_start_epoch']);assert remaining>0,'Registered 20-minute cap already elapsed'
start=perf_counter()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ball(v):return {'ball':v.str(100,more=True),'lower':v.lower().str(100,more=True),'upper':v.upper().str(100,more=True)}
def generic_dd(z):
 pos=[j for j,v in enumerate(z)if v>0]
 if len(pos)in(0,6):return fmpq(0)
 chosen,factor=(pos,2)if len(pos)<=3 else([j for j,v in enumerate(z)if v<0],-2)
 val=fmpq(0)
 for j in chosen:
  d=fmpq(1)
  for k in range(6):
   if k!=j:d*=z[j]-z[k]
  val+=fmpq(factor)/d
 return val

def laurent_dd(nodes):
 # Components: coefficient of pi/2, then coefficients of eps^-1, eps^-3, eps^-5.
 z=sorted(nodes);table=[[fmpq((v>0)-(v<0)),fmpq(0),fmpq(0),fmpq(0)]for v in z]
 for k in range(1,6):
  new=[]
  for j in range(k,6):
   if z[j]==z[j-k]:
    row=[fmpq(0)]*4
    if z[j]==0 and k%2:row[(k+1)//2]=fmpq((-1)**((k-1)//2),k)
   else:row=[(table[j-k+1][l]-table[j-k][l])/(z[j]-z[j-k])for l in range(4)]
   new.append(row)
  table=new
 return table[0]

def anchors():
 # Independent confluent polynomial interpolation, not the same DD recurrence.
 import sympy as s
 checks=[]
 for ns in [[-2,-1,-1,1,2,3],[-2,-1,0,0,1,2],[0,0,0,0,0,0]]:
  z=list(map(fmpq,ns));got=laurent_dd(z);seen={};mat=[];rhs=[[],[],[],[]]
  for x in sorted(ns):
   k=seen.get(x,0);seen[x]=k+1
   mat.append([0 if j<k else math.factorial(j)//math.factorial(j-k)*s.Integer(x)**(j-k)for j in range(6)])
   rhs[0].append((x>0)-(x<0)if k==0 else 0)
   for ii,m in enumerate([1,3,5],1):rhs[ii].append((-1)**((m-1)//2)*math.factorial(m-1)if x==0 and k==m else 0)
  M=s.Matrix(mat);expected=[M.inv()*s.Matrix(r)for r in rhs]
  assert all(str(got[l])==str(expected[l][5])for l in range(4))
  checks.append({'nodes':ns,'Laurent_coefficients':list(map(str,got)),'Hermite_interpolation_agrees':True})
 return checks

def parse_fan():
 it=iter(FAN.read_text().split());out=[]
 for _ in range(int(next(it))):
  det=int(next(it));alpha=fmpq(det);rs=[]
  for _ in range(6):
   q=int(next(it));ray=tuple(int(next(it))for _ in range(6));rs.append((q,ray));alpha/=q
  out.append((alpha,rs))
 return out

def exact_green(x,den,cones):
 tic=perf_counter();total=fmpq(0);div=[fmpq(0)]*3;terms=0;zero=0;repeat=0;nonzero=0
 for perm in itertools.permutations(range(6)):
  for mask in range(64):
   y=[(-1 if mask>>j&1 else 1)*x[perm[j]]for j in range(6)]
   for alpha,rs in cones:
    ns=sorted(fmpq(sum(y[j]*ray[j]for j in range(6)),den*q)for q,ray in rs)
    zz=any(v==0 for v in ns);rr=any(ns[j]==ns[j-1]for j in range(1,6))
    if zz or rr:
     a=laurent_dd(ns);value=a[0]
     for i in range(3):div[i]+=alpha*a[i+1]
    else:value=generic_dd(ns)
    total+=alpha*value;terms+=1;zero+=int(zz);repeat+=int(rr);nonzero+=int(value!=0)
 assert terms==967680
 assert all(v==0 for v in div),'Uncancelled divergent Laurent coefficient: singular point, no finite certificate'
 R=total/128;G=arb(R)/arb.pi()**5;assert G>0
 return {'G_times_pi5_numerator':str(R.numerator),'G_times_pi5_denominator':str(R.denominator),'G_arb384':ball(G),'complete_group_terms':terms,'nonzero_constant_terms':nonzero,'zero_node_terms':zero,'repeated_node_terms':repeat,'summed_divergent_Laurent_coefficients':list(map(str,div)),'elapsed_seconds':perf_counter()-tic}

spec=importlib.util.spec_from_file_location('b010_H_kernel',SOURCE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.FAN=FAN
x=point['candidate_numerators'];den=point['denominator'];xx=[F(v,den)for v in x]
assert xx[0]==1 and all(xx[i]>xx[i+1]>0 for i in range(5))
Ts=[F(point['T'])]+[T for T in [F(4,125),F(8,125)]if T!=F(point['T'])]
inputs=[PFILE,HERE/'gradient_results.json',HERE/'selection_results.json',HERE/'registration.md',SOURCE,FAN,Path(__file__)]
out={'status':'RUNNING_EXACT_POINT','point_index':index,'pid':os.getpid(),'point_numerators':x,'denominator':den,'preferred_T':point['T'],'new_full_passes_this_point':0,'H_rows':[],'batch_start_epoch':gr['batch_start_epoch'],'global_minimum_claimed':False,'new_C6_bound':False,'python':sys.version,'python_flint':flint.__version__,'arb_bits':384,'inputs':[{'path':str(p.resolve()),'sha256':sha(p)}for p in inputs]}
def write():
 out['elapsed_seconds']=perf_counter()-start;out['batch_elapsed_seconds']=time()-gr['batch_start_epoch']
 (HERE/f'point{index}_results.json').write_text(json.dumps(out,indent=2)+'\n')
def timeout(a,b):raise TimeoutError('Registered total batch deadline')
signal.signal(signal.SIGALRM,timeout);signal.alarm(max(1,int(remaining)))
try:
 out['scalar_confluent_anchors']=anchors();write()
 out['G']=exact_green(x,den,parse_fan());out['new_full_passes_this_point']+=1;write()
 print(json.dumps({'point':index,'G':out['G']['G_arb384']['ball'],'G_seconds':out['G']['elapsed_seconds'],'zero_terms':out['G']['zero_node_terms'],'repeated_terms':out['G']['repeated_node_terms']}),flush=True)
 G=arb(out['G']['G_arb384']['ball']);hcones=m.fan()
 for T in Ts:
  hr=m.evaluate(x,den,arb(T.numerator)/T.denominator,hcones,start,gradient=False)
  H=arb(hr['H_arb384']);K=G-H;assert K>0
  threshold=arb(T.numerator)/(T.denominator*4096);m0=arb(T.numerator)/(T.denominator*4000)
  row={'T':str(T),'H_arb384':hr['H_arb384'],'K_T':ball(K),'formal_point_cost':ball(arb(T.numerator)/(T.denominator*64*K)),
  'refutes_m_T_ge_T_over_4000':bool(K<m0),'refutes_strict_64_target':bool(K<=threshold),
  'complete_group_terms':hr['complete_group_terms'],'computed_sign_paired_terms':hr['computed_sign_paired_terms'],'elapsed_seconds':hr['elapsed_seconds']}
  out['H_rows'].append(row);out['new_full_passes_this_point']+=1;write();print(json.dumps(row),flush=True)
  if row['refutes_m_T_ge_T_over_4000']:
   out['status']='STOP_FIRST_TRUE_TARGET_COUNTEREXAMPLE';out['stop_T']=str(T);break
 else:out['status']='COMPLETE_POINT_NO_TARGET_REFUTATION'
 out['exit_code']=0;write();print(json.dumps({'status':out['status'],'pid':out['pid'],'elapsed_seconds':out['elapsed_seconds'],'batch_elapsed_seconds':out['batch_elapsed_seconds']}),flush=True)
except BaseException as e:out.update({'status':'ERROR_STOPPED','error':repr(e),'exit_code':1});write();raise
finally:signal.alarm(0)
