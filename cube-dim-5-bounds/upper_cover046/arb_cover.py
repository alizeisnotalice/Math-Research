#!/usr/bin/env python3
"""Reusable direct Arb first-jet kernel and bounded positive-mixture cover pilot."""
from __future__ import annotations
import argparse
from collections import deque
from fractions import Fraction as F
from itertools import product
from math import comb,lcm
from pathlib import Path
import gzip,hashlib,json,os,resource,sys,time
from flint import arb,fmpq,ctx
import flint

if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canonical(x)).hexdigest()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ab(x):
    x=F(x);return arb(fmpq(x.numerator,x.denominator))
def endpoint(x,upper=False):return F(str((x.upper() if upper else x.lower()).fmpq()))
def pack(x):
    if not x.is_finite():return str(x)
    return [str(endpoint(x)),str(endpoint(x,True))]
def unpack(x):
    lo,hi=map(F,x);return arb(fmpq((lo+hi).numerator,2*(lo+hi).denominator),fmpq((hi-lo).numerator,2*(hi-lo).denominator))
def atomic(p,x):
    q=p.with_suffix(p.suffix+'.tmp')
    with q.open('wb') as f:f.write(canonical(x)+b'\n');f.flush();os.fsync(f.fileno())
    os.replace(q,p)
def pointkey(x):return ','.join(map(str,x))
def numeric_quality(values):
    if not all(v.is_finite() for v in values):return False
    for v in values:
        lo,hi=map(F,pack(v))
        if hi-lo>max(abs(lo),abs(hi),F(1,10**6))/2**40:return False
    return values[0]>0


class FirstJetKernel:
    """H_epsilon and four derivatives. No frozen producer is imported.

    Integer projection eligibility is exact. All subsequent real arithmetic is
    outward Arb; partial sums are saved using exact dyadic interval endpoints.
    """
    def __init__(self,table,epsilon=F(1,64)):
        self.vertices=table['vertices'];self.simplices=table['simplices'];self.epsilon=F(epsilon)
        assert self.epsilon>0 and table['denominator']==144 and len(self.vertices)==280 and len(self.simplices)==3504

    def projections(self,point):
        point=tuple(map(F,point));assert len(point)==5
        D=lcm(*(p.denominator for p in point));n=[int(D*p) for p in point]
        k=[sum(a*b for a,b in zip(n,z)) for z in self.vertices]
        return D,k

    def eligible(self,point):
        D,k=self.projections(point)
        return {'eligible':0 not in k and len(set(k))==280,'zero_count':k.count(0),
                'distinct_count':len(set(k)),'common_denominator':str(D)}

    def evaluate(self,point,bits=256,partial=None,deadline=None,stride=128):
        """Return {complete, values/partial, ...}; a deadline never omits a cone."""
        began=time.monotonic();ctx.prec=bits;point=tuple(map(F,point));D,k=self.projections(point)
        assert 0 not in k and len(set(k))==280,'generic-center-only DD API'
        scaled=144*D*self.epsilon;A,B=scaled.numerator,scaled.denominator
        slopes=[tuple(D*z for z in v[1:]) for v in self.vertices]
        jets=[]
        for node,l in zip(k,slopes):
            w=B*B*node*node+A*A
            jets.append([arb(w).log()/2,*[arb(s*B*B*node)/w for s in l]])
        total=[arb(0) for _ in range(5)];absolute=arb(0);start=0
        if partial:
            assert partial['point']==list(map(str,point)) and partial['bits']==bits
            start=partial['next_cone'];total=list(map(unpack,partial['total']));absolute=unpack(partial['absolute'])
        for ci in range(start,3504):
            s=self.simplices[ci];ids=sorted(s['vertices'],key=lambda i:k[i]);table=[]
            for right,idx in enumerate(ids):
                row=[jets[idx]]
                for order in range(1,right+1):
                    left=ids[right-order];d=k[idx]-k[left]
                    delta=[a-b for a,b in zip(slopes[idx],slopes[left])]
                    P=row[-1];Q=table[right-1][order-1]
                    r0=(P[0]-Q[0])/d
                    row.append([r0,*[(P[j+1]-Q[j+1]-delta[j]*r0)/d for j in range(4)]])
                table.append(row)
            jet=table[-1][-1];det=s['det_numerator']
            for j in range(5):total[j]-=det*jet[j]
            absolute+=abs(det*jet[0]);done=ci+1
            if deadline is not None and done%stride==0 and time.monotonic()>=deadline:
                return {'complete':False,'partial':{'point':list(map(str,point)),'bits':bits,'next_cone':done,
                  'total':list(map(pack,total)),'absolute':pack(absolute)},'seconds':time.monotonic()-began}
        factor=arb(D)**4/(144*(2*arb.pi())**5);values=[factor*v for v in total]
        return {'complete':True,'values':values,'seconds':time.monotonic()-began,
          'cancellation_ratio':str(absolute*factor/abs(values[0])) if not values[0].contains(0) else None,
          'minimum_integer_projection_gap':str(min(b-a for a,b in zip(sorted(k),sorted(k)[1:]))),
          'common_denominator':str(D),'all_cones':3504,'bits':bits}


def roots():
    return [[i,j,k,l] for i in range(32) for j in range(i+1) for k in range(j+1) for l in range(k+1)]
def nodekey(node):return str(node[0])+':'+node[1]
def box_of(node,rootlist):
    lo=[F(i,32) for i in rootlist[node[0]]];hi=[x+F(1,32) for x in lo]
    for bit in node[1]:
        side=[b-a for a,b in zip(lo,hi)];axis=max(range(4),key=lambda j:side[j]);mid=(lo[axis]+hi[axis])/2
        if bit=='0':hi[axis]=mid
        else:assert bit=='1';lo[axis]=mid
    return lo,hi
def choose_anchor(box,kernel):
    lo,hi=box;mid=[(a+b)/2 for a,b in zip(lo,hi)];eta=min(b-a for a,b in zip(lo,hi))/4096
    candidates=[mid]+[[x+eta*beta**j for j,x in enumerate(mid)] for beta in (F(1,8),F(1,16),F(1,32))]
    diagnostics=[]
    for i,tail in enumerate(candidates):
        assert all(a<=x<=b for a,x,b in zip(lo,tail,hi))
        point=(F(1),*tail);test=kernel.eligible(point);diagnostics.append({'candidate':i,**test})
        if test['eligible']:return point,i,diagnostics
    return None,None,diagnostics
def certificate(box,record,epsilon,target):
    Hlo=F(record['values'][0][0]);g=[unpack(x)/unpack(record['values'][0]) for x in record['values'][1:]]
    gd=[tuple(map(F,pack(x))) for x in g];Delta=endpoint((ab(Hlo)/ab(target)).log())
    point=list(map(F,record['point']));lo,hi=box;rows=[]
    for bits in product((0,1),repeat=4):
        u=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)]
        b=sum((min(a*v,b*v) for (a,b),v in zip(gd,u)),F())
        q=sum((v*v for v in u),F())/4-sum(u,F())**2/36
        margin=Delta+b-q/(epsilon*epsilon);rows.append((margin,bits))
    margin,bits=min(rows)
    return {'passed':margin>0,'minimum_log_margin':str(margin),'worst_corner_bits':list(bits),
            'Delta_lower':str(Delta),'loggradient_directed':[[str(a),str(b)] for a,b in gd],
            'all16_corners':True,'quadratic':'Q=I4/4-J4/36','anchor':pointkey(point)}
def check_tree(state):
    leaves={nodekey([i,'']) for i in range(len(state['roots']))}
    for split in state['splits']:
        key=nodekey(split['node']);assert key in leaves;leaves.remove(key)
        for child in split['children']:assert nodekey(child) not in leaves;leaves.add(nodekey(child))
    pending=[nodekey(n) for n in state['queue']];passed=[nodekey(n['node']) for n in state['passed']]
    assert len(pending)==len(set(pending)) and len(passed)==len(set(passed))
    assert set(pending).isdisjoint(passed) and leaves==set(pending)|set(passed)
    assert len(pending)+len(passed)==len(state['roots'])+len(state['splits'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream-root',type=Path,required=True)
    parser.add_argument('--reference-root',type=Path,default=Path(__file__).resolve().parent.parent/'upper_harnack043')
    parser.add_argument('--work-dir',type=Path,required=True);parser.add_argument('--resume',action='store_true')
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args();began=time.monotonic();base=Path(__file__).resolve().parent
    reg=json.loads((base/'registration.json').read_text());epsilon=F(reg['epsilon']);target=F(reg['target'])
    for rel,h in reg['upstream_read_only_inputs'].items():assert sha(args.upstream_root/rel)==h
    for rel,h in reg['reference043_read_only_inputs'].items():assert sha(args.reference_root/rel)==h
    sources={'registration.json':sha(base/'registration.json'),'arb_cover.py':sha(Path(__file__))}
    fingerprint=digest({'sources':sources,'inputs':reg['upstream_read_only_inputs'],'reference':reg['reference043_read_only_inputs']})
    work=args.work_dir.resolve();assert not any((p/'.git').exists() for p in (work,*work.parents));work.mkdir(parents=True,exist_ok=True)
    checkpoint=work/'state.json'
    if checkpoint.exists():
        assert args.resume,'existing state requires --resume';state=json.loads(checkpoint.read_text());assert state['fingerprint']==fingerprint
        assert state['status'] not in ('fixed_epsilon_H_counterexample','complete_cover','qualification_failed','numerical_diagnostic','no_generic_anchor')
    else:
        assert not args.resume;initial=roots();assert len(initial)==comb(35,4)==52360
        queue=[[i,''] for i in range(len(initial))]
        state={'fingerprint':fingerprint,'status':'initialized','roots':initial,'queue':queue,'splits':[],'passed':[],
          'initial_root_sha256':digest(initial),'initial_queue_sha256':digest(queue),'cache':{},'qualification':None,
          'active':None,'segments':[],'total_new_anchor_starts':0,'total_kernel_attempts_completed':0,'cache_hits':0,
          'numerical_diagnostics':[],'geometry_candidates_rejected':0,'source_hashes':sources}
    def save(status):
        state['status']=status;check_tree(state);atomic(checkpoint,state)
    save(state['status'])
    table=json.loads(gzip.decompress((args.upstream_root/'astra_vertex017/boundary_triangulation.json.gz').read_bytes()))
    kernel=FirstJetKernel(table,epsilon);limit=reg['finite_batch_limits'];deadline=began+limit['seconds'];newstarts=0;newcompleted=0;calls=0;status='partial_cover'
    reference=json.loads((args.reference_root/'results.json').read_text())
    # Exactly one old-coordinate qualification object; precision retry is not a new point.
    if state['qualification'] is None:
        qpoint=tuple(map(F,reg['qualification']['point']));qsteps=[]
        if state['active'] is not None:
            assert state['active']['kind']=='qualification';qsteps=state['active']['attempts']
        for pi,bits in enumerate(reg['qualification']['precision_ladder_bits']):
            if pi<len(qsteps):continue
            active=state['active'] or {'kind':'qualification','attempts':qsteps,'partial':None}
            state['active']=active
            got=kernel.evaluate(qpoint,bits,active['partial'],deadline,limit['kernel_checkpoint_stride_cones']);calls+=1
            if not got['complete']:
                active['partial']=got['partial'];status='time_checkpoint';break
            old=[arb(reference['H_center_ball']),*map(arb,reference['four_tail_gradient_H_balls'])]
            checks=[v.contains(w) for v,w in zip(got['values'],old)];qualified=all(checks) and numeric_quality(got['values'])
            qsteps.append({'bits':bits,'seconds':got['seconds'],'values':list(map(pack,got['values'])),
                           'contains_five043_balls':checks,'numeric_quality':numeric_quality(got['values'])})
            active['partial']=None;active['attempts']=qsteps
            if qualified:
                state['qualification']={'passed':True,'point':list(map(str,qpoint)),'attempts':qsteps};state['active']=None;save('qualified');break
        if state['qualification'] is None and status!='time_checkpoint':status='qualification_failed'
    queue=deque(state['queue'])
    while state['qualification'] is not None and queue and status=='partial_cover':
        if time.monotonic()>=deadline:status='time_checkpoint';break
        if state['active'] is None:
            if newstarts>=limit['max_new_anchor_starts']:status='anchor_checkpoint';break
            node=queue[0];box=box_of(node,state['roots']);point,which,tests=choose_anchor(box,kernel)
            state['geometry_candidates_rejected']+=sum(not t['eligible'] for t in tests)
            if point is None:
                state['numerical_diagnostics'].append({'node':node,'geometry_tests':tests});status='no_generic_anchor';break
            key=pointkey(point)
            if key in state['cache']:record=state['cache'][key];state['cache_hits']+=1
            else:
                state['active']={'kind':'pilot','node':node,'point':list(map(str,point)),
                  'candidate':which,'geometry_tests':tests,'precision_index':0,'partial':None,'attempts':[]}
                state['total_new_anchor_starts']+=1;newstarts+=1;record=None
        else:
            active=state['active'];assert active['kind']=='pilot' and active['node']==queue[0]
            node=active['node'];box=box_of(node,state['roots']);record=None
        if record is None:
            active=state['active'];pi=active['precision_index'];bits=reg['pilot_precision_ladder_bits'][pi]
            got=kernel.evaluate(active['point'],bits,active['partial'],deadline,limit['kernel_checkpoint_stride_cones']);calls+=1
            if not got['complete']:
                active['partial']=got['partial'];status='time_checkpoint';break
            values=got['values'];packed=list(map(pack,values));quality=numeric_quality(values)
            state['total_kernel_attempts_completed']+=1
            active['attempts'].append({'bits':bits,'seconds':got['seconds'],'values':packed,'numeric_quality':quality,
                  'cancellation_ratio':got['cancellation_ratio'],'common_denominator':got['common_denominator'],
                  'minimum_integer_projection_gap':got['minimum_integer_projection_gap']})
            active['partial']=None
            if values[0].is_finite() and endpoint(values[0],True)<target:
                state['counterexample']={'point':active['point'],'node':node,'H_interval':packed[0],
                   'strict_H_upper_below_target':True,'target':str(target),'bits':bits,
                   'scope':'fixed epsilon H target only; no full-G or C5 failure'}
                status='fixed_epsilon_H_counterexample';newcompleted+=1;break
            separated=values[0].is_finite() and endpoint(values[0])>target
            if not quality or not separated:
                if pi+1<len(reg['pilot_precision_ladder_bits']):active['precision_index']+=1;save('precision_retry');continue
                state['numerical_diagnostics'].append({'node':node,'point':active['point'],'attempts':active['attempts']});status='numerical_diagnostic';break
            record={'point':active['point'],'values':packed,'bits':bits,'candidate':active['candidate'],
                    'geometry_tests':active['geometry_tests'],'attempts':active['attempts']}
            state['cache'][pointkey(active['point'])]=record;state['active']=None;newcompleted+=1
        cert=certificate(box,record,epsilon,target);head=queue.popleft();assert head==node
        if cert['passed']:state['passed'].append({'node':node,'certificate':cert})
        else:
            lo,hi=box;axis=max(range(4),key=lambda j:hi[j]-lo[j]);children=[[node[0],node[1]+'0'],[node[0],node[1]+'1']]
            state['splits'].append({'node':node,'axis':axis,'children':children,'failed_certificate':cert});queue.extend(children)
        state['queue']=list(queue);save('partial_cover')
        print(f"anchors {newcompleted}/{newstarts}; pass {len(state['passed'])}; split {len(state['splits'])}; queue {len(queue)}; {time.monotonic()-began:.3f}s",flush=True)
    if not queue and status=='partial_cover':status='complete_cover'
    state['queue']=list(queue);state['segments'].append({'pid':os.getpid(),'seconds':time.monotonic()-began,
      'new_anchor_starts':newstarts,'new_anchors_completed':newcompleted,'kernel_calls':calls,'stop':status})
    save(status)
    certs=[F(x['certificate']['minimum_log_margin']) for x in state['passed']]
    tested=certs+[F(x['failed_certificate']['minimum_log_margin']) for x in state['splits']]
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result={'status':status,'model':'alpha1','epsilon':str(epsilon),'target':str(target),'global_g_or_C5_claim':False,
      'qualification':state['qualification'],'initial_roots':len(state['roots']),'initial_root_sha256':state['initial_root_sha256'],
      'initial_queue_sha256':state['initial_queue_sha256'],'total_new_anchor_starts':state['total_new_anchor_starts'],
      'completed_cached_anchors':len(state['cache']),'total_kernel_attempts_completed':state['total_kernel_attempts_completed'],
      'passed_leaves':len(state['passed']),'split_nodes':len(state['splits']),'remaining_queue_leaves':len(queue),
      'initial_roots_never_processed':len(state['roots'])-len({x['node'][0] for x in state['splits']+state['passed']}),
      'geometry_candidates_rejected':state['geometry_candidates_rejected'],'cache_hits':state['cache_hits'],
      'worst_passed_log_margin':str(min(certs)) if certs else None,'worst_tested_log_margin':str(min(tested)) if tested else None,
      'cover_tree_sha256':digest({'roots':state['roots'],'splits':state['splits'],'passed':state['passed'],'queue':state['queue']}),
      'remaining_queue_sha256':digest(state['queue']),'state':{'path':str(checkpoint),'sha256':sha(checkpoint),'bytes':checkpoint.stat().st_size,'fingerprint':fingerprint},
      'counterexample':state.get('counterexample'),'diagnostics':state['numerical_diagnostics'],
      'runtime':{'wall_seconds':time.monotonic()-began,'peak_rss_bytes':rss if sys.platform=='darwin' else rss*1024,
          'flint_version':flint.__version__,'single_process':True,'segments':state['segments']},
      'source_hashes':sources,'upstream_hashes':reg['upstream_read_only_inputs'],'reference043_hashes':reg['reference043_read_only_inputs']}
    args.output_dir.mkdir(parents=True,exist_ok=True);out=args.output_dir/'results.json';atomic(out,result)
    print(json.dumps({key:result[key] for key in ('status','total_new_anchor_starts','passed_leaves','split_nodes','remaining_queue_leaves','runtime','state')},indent=2))


if __name__=='__main__':main()
