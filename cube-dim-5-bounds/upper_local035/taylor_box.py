#!/usr/bin/env python3
"""One full-time Cauchy Green Taylor box, exact rational log stencils.

No frozen017/018 code is executed. The3504-simplex table and stored center
ball are read-only, SHA-locked inputs. All resumable stencils stay outside Git.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from itertools import product
from math import factorial, comb, lcm, prod
from pathlib import Path
import gzip, hashlib, json, os, resource, sys, tempfile, time
import flint
from flint import arb, fmpq, ctx

if hasattr(sys,'set_int_max_str_digits'): sys.set_int_max_str_digits(0)


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rq(s):
    s=str(s)
    if '/' in s:
        a,b=s.split('/');return fmpq(int(a),int(b))
    return fmpq(int(s))
def aq(q): return arb(q)
def endpoint(x,upper=False): return str((x.upper() if upper else x.lower()).fmpq())
def atomic(path,data):
    temporary=path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w') as f:
        json.dump(data,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(temporary,path)


def atanh_log_upper(u,N=16):
    y=(u-1)/(u+1);assert 0<=y<=fmpq(1,3)
    total=fmpq(0);power=y
    for j in range(N): total+=2*power/(2*j+1);power*=y*y
    return total+2*power/((2*N+1)*(1-y*y))


LOG2_UPPER=atanh_log_upper(fmpq(2))
def log_ratio_upper(a,b):
    # log(b/a), integers0<a<=b. No floating logarithm or fitting.
    assert 0<a<=b
    e=max(0,b.bit_length()-a.bit_length())
    if a*(1<<e)>b:e-=1
    u=fmpq(b,a*(1<<e));assert 1<=u<2
    return e*LOG2_UPPER+atanh_log_upper(u)


def time_absolute_upper(nodes):
    a=sorted(map(abs,nodes));assert a[0]>0
    answer=fmpq(a[0],prod(a))
    answer+=log_ratio_upper(a[0],a[1])/prod(a[1:])
    for j in range(2,5):
        answer+=(fmpq(1,a[j-1]**(j-1))-fmpq(1,a[j]**(j-1)))/((j-1)*prod(a[j:]))
    answer+=fmpq(1,4*a[4]**4)
    return answer


def simplex_stencil(ids,vertices,k,D,alphas,parents,alpha_factor):
    ids=sorted(ids,key=lambda i:k[i]);nodes=[k[i] for i in ids]
    assert len(set(nodes))==5 and 0 not in nodes
    slopes=[tuple(D*z for z in vertices[i][1:]) for i in ids]
    table=[];nalpha=len(alphas)
    for right,node in enumerate(nodes):
        jet=[[fmpq(0) for _ in range(6)] for _ in alphas]
        jet[0][right+1]=fmpq(1)
        for ai,alpha in enumerate(alphas[1:],1):
            value=alpha_factor[ai]
            for slope,power in zip(slopes[right],alpha):
                if power:value*=fmpq(slope,node)**power
            jet[ai][0]=value
        row=[jet]
        for order in range(1,right+1):
            left=right-order;d=node-nodes[left]
            dl=[a-b for a,b in zip(slopes[right],slopes[left])]
            P=row[-1];Q=table[right-1][order-1]
            quotient=[]
            for ai in range(nalpha):
                coefficient=[]
                for symbol in range(6):
                    numerator=P[ai][symbol]-Q[ai][symbol]
                    for j,previous in parents[ai]:
                        numerator-=dl[j]*quotient[previous][symbol]
                    coefficient.append(numerator/d)
                quotient.append(coefficient)
            row.append(quotient)
        table.append(row)
    return table[-1][-1],[abs(n) for n in nodes],nodes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream-root',type=Path,required=True,help='Read-only project directory containing astra_vertex017/ and astra_green018/')
    parser.add_argument('--work-dir',type=Path)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args();began=time.monotonic();root=Path(__file__).resolve().parent
    registration_path=root/'registration.json';reg=json.loads(registration_path.read_text())
    ctx.prec=reg['precision_bits'];M=reg['Taylor_total_degree'];h=rq(reg['one_closed_box']['four_tail_halfwidth'])
    for relative,digest in reg['read_only_inputs'].items():
        assert sha(args.upstream_root/relative)==digest,(relative,'input hash mismatch')
    source_hash=sha(Path(__file__));reg_hash=sha(registration_path)
    fingerprint=hashlib.sha256(json.dumps({'source':source_hash,'registration':reg_hash,
        'inputs':reg['read_only_inputs']},sort_keys=True).encode()).hexdigest()
    work=args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix='upper-local035-'))
    work.mkdir(parents=True,exist_ok=True)
    if any((p/'.git').exists() for p in (work,*work.parents)):raise ValueError('work/cache must be outside Git')
    checkpoint=work/'checkpoint.json';state_path=work/'state.json'
    if checkpoint.exists() and not args.resume:raise ValueError('existing checkpoint requires --resume')
    table=json.loads(gzip.decompress((args.upstream_root/'astra_vertex017/boundary_triangulation.json.gz').read_bytes()))
    assert table['denominator']==144 and len(table['vertices'])==280 and len(table['simplices'])==3504
    vertices=table['vertices'];simplices=table['simplices']
    point=tuple(map(F,reg['one_closed_box']['center']));D=lcm(*(x.denominator for x in point))
    assert D==104729 and point[0]==1
    center_n=[int(x*D) for x in point]
    k=[sum(a*b for a,b in zip(center_n,z)) for z in vertices]
    assert 0 not in k and len(set(k))==280
    slopes_norm=[D*sum(map(abs,z[1:])) for z in vertices]
    q=[h*fmpq(L,abs(n)) for L,n in zip(slopes_norm,k)]
    assert max(q)<1,'fixed box fails nonzero projection guard; no adjustment permitted'
    alphas=sorted((a for a in product(range(M+1),repeat=4) if sum(a)<=M),key=lambda a:(sum(a),a))
    assert len(alphas)==35 and alphas[0]==(0,0,0,0)
    index={a:i for i,a in enumerate(alphas)};parents=[];alpha_factor=[]
    for a in alphas:
        parents.append([(j,index[tuple(t-int(j==l) for l,t in enumerate(a))]) for j in range(4) if a[j]])
        n=sum(a);alpha_factor.append(fmpq(0) if n==0 else fmpq((-1)**(n-1)*factorial(n-1),prod(factorial(v) for v in a)))
    stencil=[{} for _ in alphas];remainder=arb(0);next_cone=0;segments=[]
    if args.resume:
        saved=json.loads(checkpoint.read_text());assert saved['fingerprint']==fingerprint
        next_cone=saved['next_simplex'];segments=saved['segments']
        stencil=[{int(key):rq(value) for key,value in row.items()} for row in saved['stencil']]
        remainder=aq(rq(saved['remainder_sum_upper']))
    def save(next_index,status):
        elapsed=time.monotonic()-began
        payload={'fingerprint':fingerprint,'status':status,'next_simplex':next_index,
          'stencil':[{str(key):str(value) for key,value in row.items()} for row in stencil],
          'remainder_sum_upper':endpoint(remainder,True),
          'segments':segments+[{'pid':os.getpid(),'seconds':elapsed,'end_simplex':next_index}]}
        atomic(checkpoint,payload)
        atomic(state_path,{'fingerprint':fingerprint,'status':status,'next_simplex':next_index,
                           'PID':os.getpid(),'seconds_this_process':elapsed,'checkpoint':str(checkpoint)})
    for ci in range(next_cone,len(simplices)):
        simplex=simplices[ci];det=int(simplex['det_numerator']);ids=simplex['vertices']
        jet,log_keys,nodes=simplex_stencil(ids,vertices,k,D,alphas,parents,alpha_factor)
        for ai,row in enumerate(jet):
            global_row=stencil[ai]
            for symbol,value in enumerate(row):
                if value==0:continue
                key=0 if symbol==0 else log_keys[symbol-1]
                combined=global_row.get(key,fmpq(0))-det*value
                if combined==0:global_row.pop(key,None)
                else:global_row[key]=combined
        qc=max(q[i] for i in ids)
        tail=(1-qc)**(-5)-sum((fmpq(comb(n+4,4))*qc**n for n in range(M+1)),fmpq(0))
        assert tail>=0
        remainder+=aq(det*time_absolute_upper(nodes)*tail)
        done=ci+1
        if done%reg['checkpoint_stride_full_simplices']==0 or done==len(simplices):
            save(done,'complete_stencils' if done==3504 else 'checkpoint')
            print(f'completed {done}/3504 simplices; elapsed {time.monotonic()-began:.3f}s',flush=True)
        if time.monotonic()-began>=reg['checkpoint_seconds'] and done<3504:
            save(done,'checkpoint');raise SystemExit(75)
    # Evaluate the globally collected EXACT rational linear log stencils.
    logkeys=sorted({key for row in stencil for key in row if key})
    logs={key:arb(key).log() for key in logkeys}
    factor=arb(D)**4/(144*(2*arb.pi())**5);coefficients=[]
    for row in stencil:
        value=arb(0)
        for key,c in row.items():value+=aq(c)*(logs[key] if key else 1)
        coefficients.append(factor*value)
    old=json.loads((args.upstream_root/'astra_green018/batch_result.json').read_text())['selected_value']
    assert tuple(map(F,old['point']))==point
    old_ball=arb(old['G_ball']);constant=coefficients[0]
    assert not (constant.upper()<old_ball.lower() or old_ball.upper()<constant.lower()),'center compatibility failed'
    variation=arb(0)
    for a,c in zip(alphas[1:],coefficients[1:]):variation+=abs(c)*aq(h**sum(a))
    full_remainder=factor*aq(rq(endpoint(remainder,True)))
    lower=constant-variation-full_remainder;upper=constant+variation+full_remainder
    target=aq(rq(reg['one_closed_box']['target_lower_bound']))
    passed=lower>target
    final_stencil_size=sum(len(row) for row in stencil)
    bit_max=0
    for row in stencil:
        for c in row.values():
            parts=str(c).split('/')
            bit_max=max(bit_max,*(abs(int(s)).bit_length() for s in parts))
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result={'status':'certified_one_closed_spatial_box' if passed else 'fixed_box_bound_not_separated',
      'model':'alpha=1, full3504-simplex Fourier Green, full time0..infinity',
      'global_g_or_C5_claim':False,'box':reg['one_closed_box'],'Taylor_total_degree':M,
      'precision_bits':ctx.prec,'all_simplices_processed':3504,'all_vertices_qualified':280,
      'maximum_vertex_relative_projection_shift_exact':str(max(q)),
      'minimum_center_abs_integer_projection':min(map(abs,k)),
      'center_projection_distinctness_only':True,
      'repeated_nonzero_projections_inside_box':'Handled by the complete-time resolvent Taylor proof, not ordinary DD substitution.',
      'Taylor_coefficients':[{'alpha':list(a),'ball':str(c)} for a,c in zip(alphas,coefficients)],
      'center_G_ball':str(constant),'stored018_center_G_ball':old['G_ball'],'center_balls_overlap':True,
      'Taylor_polynomial_variation_upper_ball':str(variation),
      'complete_time_absolute_remainder_upper_ball':str(full_remainder),
      'G_uniform_lower_dyadic':endpoint(lower,False),'G_uniform_upper_dyadic':endpoint(upper,True),
      'G_uniform_lower_display':str(lower),'G_uniform_upper_display':str(upper),
      'target_certified':bool(passed),'exact_rational_stencil_entries':final_stencil_size,
      'maximum_exact_stencil_integer_bits':bit_max,
      'checkpoint':{'path':str(checkpoint),'sha256':sha(checkpoint),'bytes':checkpoint.stat().st_size,
                    'fingerprint':fingerprint,'segments':segments+[{'pid':os.getpid(),'seconds':time.monotonic()-began,'end_simplex':3504}]},
      'runtime':{'PID':os.getpid(),'wall_seconds_this_process':time.monotonic()-began,
        'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,'flint_version':flint.__version__,
        'boxes':1,'mesh_cells':0,'time_cutoffs':0},
      'source_sha256':{'registration.json':reg_hash,'taylor_box.py':source_hash,**reg['read_only_inputs']}}
    args.output_dir.mkdir(parents=True,exist_ok=True);destination=args.output_dir/'results.json'
    atomic(destination,result)
    print(json.dumps({'status':result['status'],'output':str(destination),'sha256':sha(destination),
      'G_lower':str(lower),'variation':str(variation),'remainder':str(full_remainder),
      'maximum_q':str(max(q)),'runtime':result['runtime']},indent=2))


if __name__=='__main__':main()
