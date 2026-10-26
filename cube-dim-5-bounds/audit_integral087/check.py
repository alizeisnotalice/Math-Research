#!/usr/bin/env python3
"""087 independent audit: standard-library exact arithmetic only; no H/Arb calls."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement,permutations,product
from math import isqrt,lcm,prod,factorial,ceil
import argparse,ast,gzip,hashlib,json,os,sys,time

HERE=Path(__file__).resolve().parent
N=1<<160
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def read(p):
    p=Path(p)
    return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())


def atomic(p,o):
    t=p.with_suffix('.tmp');t.write_text(json.dumps(o,sort_keys=True,separators=(',',':'))+'\n');os.replace(t,p)


def manifest_check(p):
    m=read(p);files=m['files']
    if isinstance(files,dict):items=[dict(v,name=k) for k,v in files.items()]
    else:items=files
    for row in items:
        f=Path(p).parent/row.get('name',row.get('path'))
        assert sha(f)==row['sha256'],str(f)
    return len(items)


def orbit(i):return set(permutations(i))
def weight(i):return len(orbit(i))
def key(level,i):return str(level)+':'+','.join(map(str,i))
def idx(k):return tuple(map(int,k.split(':')[1].split(',')))
def qform(u):return sum(x*x for x in u)/4-sum(u)**2/36


def root4(n):
    a=isqrt(isqrt(n));assert a**4<=n<(a+1)**4;return a


def chord_check(b,level):
    li,ui,mi=[b[k] for k in ('L_i','U_i','M_i')]
    assert all(isinstance(x,int) for x in (li,ui,mi,b['integral_lower_i'],b['integral_upper_i']))
    assert 0<li<=mi<=ui<=3*N and li<ui and b['positive_bands']==73
    a=root4(li**5);z=root4(ui**5)
    # (li/2^160)^(5/4) lies in [a,a+1]/2^200. No transcendental library.
    shift=1<<(40+4*level)
    assert b['integral_lower_i']*shift<=a
    uppernum=(ui-mi)*(a+1)+(mi-li)*(z+1)
    assert b['integral_upper_i']*(ui-li)*shift>=uppernum
    assert 0<=b['integral_lower_i']<=b['integral_upper_i']


def geometry(point,level,indices):
    tail=list(map(F,point))[1:];width=F(1,1<<level)
    lo=[i*width-c for i,c in zip(indices,tail)]
    hi=[(i+1)*width-c for i,c in zip(indices,tail)]
    assert all(l<0<u for l,u in zip(lo,hi))
    mu=[(l+u)/2 for l,u in zip(lo,hi)]
    second=[(l*l+l*u+u*u)/3 for l,u in zip(lo,hi)]
    eta=sum(second)/4-(sum(second)+sum(mu[i]*mu[j] for i in range(4) for j in range(4) if i!=j))/36
    return {'lo':list(map(str,lo)),'hi':list(map(str,hi)),'mu':list(map(str,mu)),
            'eta':str(eta),'B2':str(max(qform(v) for v in product(*zip(lo,hi)))),
            'h_inf':str(max(map(abs,lo+hi))),'lambda':str(max(map(abs,mu))),
            'volume':str(width**4)}


def projections(point,vertices):
    point=list(map(F,point));assert point[0]==1
    D=lcm(*(x.denominator for x in point));a=[int(x*D) for x in point]
    ns=[sum(x*z for x,z in zip(a,v)) for v in vertices]
    assert len(set(ns))==280
    return D,ns


def partition(parent,child_indices,boundary=False):
    ordered={tuple(2*a+b for a,b in zip(p,bits)) for p in orbit(parent) for bits in product((0,1),repeat=4)}
    expanded=set()
    for c in child_indices:
        assert tuple(a//2 for a in c)==tuple(parent)
        o=orbit(c);assert not expanded&o;expanded|=o
    assert expanded==ordered
    assert sum(weight(c) for c in child_indices)==16*weight(parent)
    if boundary:
        for offsets in product((0,1,2),repeat=4):
            p=tuple(2*i+s for i,s in zip(parent,offsets))
            candidates=[tuple(j for j in (x-1,x) if 2*a<=j<=2*a+1) for a,x in zip(parent,p)]
            assert all(candidates) and any(v in expanded for v in product(*candidates))
    return len(ordered)


def children(p):
    # Independent group orbit construction, rather than sorting representative bits.
    all_ordered={tuple(2*a+b for a,b in zip(o,bits)) for o in orbit(p) for bits in product((0,1),repeat=4)}
    return sorted({tuple(sorted(c,reverse=True)) for c in all_ordered})


def fifth_root_floor(n):
    x=1<<((n.bit_length()+4)//5)
    while True:
        y=(4*x+n//(x**4))//5
        if y>=x:break
        x=y
    while x**5>n:x-=1
    while (x+1)**5<=n:x+=1
    return x


def full_I_cap(UH,E,cap):
    scale=1<<256;z=UH**4
    a=fifth_root_floor(z.numerator*scale**5//z.denominator)
    B=F(a+1,scale)+E
    assert cap**4>=B**5
    return str(B)


def synthetic_cuts():
    root=(0,0,0,0);cs=children(root);grand={c:children(c) for c in cs}
    assert len(cs)==5 and sum(map(len,grand.values()))==35
    cuts=[[(0,root)]]
    for bits in product((0,1),repeat=5):
        cut=[]
        for c,split in zip(cs,bits):cut+=([(2,z) for z in grand[c]] if split else [(1,c)])
        cuts.append(cut)
    assert len(cuts)==33
    def truth(level,i):
        m=1<<level;v=F(weight(i),m**4)
        return v*(1+sum(F(a*a+a*(a+1)+(a+1)**2,3*m*m) for a in i))
    for cut in cuts:
        assert sum((F(weight(i),1<<(4*level)) for level,i in cut),F(0))==1
        assert sum((truth(level,i) for level,i in cut),F(0))==F(7,3)
    cases=[]
    allnodes=[(0,root)]+[(1,c) for c in cs]+[(2,z) for c in cs for z in grand[c]]
    for case in range(3):
        vals={}
        for level,i in allnodes:
            v=truth(level,i)
            delta=([F(1,100),F(1,50),F(1,25)][level] if case==0 else
                   [F(1,10),F(1,100),F(1,10000)][level] if case==1 else
                   F(1+(sum(i)%3),10**(1+(level+sum(i))%4)))
            vals[(level,i)]=(max(F(0),v-delta),v+delta)
        lows=[sum((vals[n][0] for n in cut),F(0)) for cut in cuts]
        ups=[sum((vals[n][1] for n in cut),F(0)) for cut in cuts]
        child_eff=[]
        for c in cs:
            lower=max(vals[(1,c)][0],sum((vals[(2,z)][0] for z in grand[c]),F(0)))
            upper=min(vals[(1,c)][1],sum((vals[(2,z)][1] for z in grand[c]),F(0)))
            child_eff.append((lower,upper))
        eff=(max(vals[(0,root)][0],sum(x[0] for x in child_eff)),min(vals[(0,root)][1],sum(x[1] for x in child_eff)))
        assert eff==(max(lows),min(ups)) and eff[0]<=F(7,3)<=eff[1]
        cases.append({'case':case,'cut_count':33,'effective_lower':str(eff[0]),'effective_upper':str(eff[1]),'truth':'7/3'})
    return cases


def main():
    start=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');args=ap.parse_args()
    reg=read(HERE/'registration.json');ext=read(HERE/'scope_extension.json')
    assert sha(HERE/'registration.json')==ext['original_registration_sha256']
    assert sha(HERE/'audit_protocol.md')==reg['protocol_sha256']==ext['original_protocol_sha256']
    assert sha(HERE/'sample082.json.gz')==reg['sample082_sha256']
    sources=reg['source_files'];extra=ext['extra_files']
    for s in list(sources.values())+list(extra.values()):assert sha(s['path'])==s['sha256'],s['path']
    path={k:Path(v['path']) for k,v in sources.items()}
    manifest_counts={k:manifest_check(v) for k,v in path.items() if k.endswith('manifest')}
    manifest_counts['082_final']=manifest_check(extra['manifest.json']['path'])
    fp=hashlib.sha256((sha(__file__)+sha(HERE/'registration.json')+sha(HERE/'scope_extension.json')).encode()).hexdigest()
    cpfile=HERE/'checkpoint.json';cp={'fingerprint':fp,'completed_roots':0,'lower_units':0,'upper_units':0,'weight':0}
    if cpfile.exists():assert args.resume;cp=read(cpfile);assert cp['fingerprint']==fp
    roots=read(path['080_roots']);compact=read(path['080_cells']);cols=compact['columns'];rows=[dict(zip(cols,r)) for r in compact['rows']]
    assert len(rows)==len(roots)==52360 and F(compact['denominator'])==N
    catalog=sorted((tuple(reversed(x)) for x in combinations_with_replacement(range(32),4)))
    assert len(catalog)==52360 and sorted(r['rootid'] for r in roots)==list(range(52360))
    vertices=read(path['017_table'])['vertices'];table=read(path['017_table'])
    assert len(vertices)==280 and len(table['simplices'])==3504 and table['denominator']==144
    detsum=sum(c['det_numerator'] for c in table['simplices']);assert all(c['det_numerator']>0 and len(set(c['vertices']))==5 for c in table['simplices'])
    previous=None
    for rank,(root,b) in enumerate(zip(roots,rows)):
        score=F(root['weight'])**4*F(root['cached_H1over64']['upper'])**5
        if rank>=cp['completed_roots']:
            assert root['rank']==b['rank']==rank and root['rootid']==b['rootid']
            i=tuple(root['indices']);assert i==catalog[root['rootid']]
            w=weight(i);assert root['weight']==b['weight']==w
            assert 0<F(root['cached_H1over64']['lower'])<=F(root['cached_H1over64']['upper'])
            geometry(root['point'],5,i);projections(root['point'],vertices);chord_check(b,5)
            assert previous is None or score<previous[0] or (score==previous[0] and root['rootid']>previous[1])
            cp['lower_units']+=w*b['integral_lower_i'];cp['upper_units']+=w*b['integral_upper_i'];cp['weight']+=w;cp['completed_roots']=rank+1
            if (rank+1)%4096==0 or rank==52359:
                cp['elapsed_seconds']=time.monotonic()-start;atomic(cpfile,cp);print('080 checked',rank+1,'seconds',round(time.monotonic()-start,3),flush=True)
        previous=(score,root['rootid'])
    summary80=read(path['080_summary'])
    assert cp['weight']==32**4 and F(cp['lower_units'],N)==F(summary80['H_integral_lower']) and F(cp['upper_units'],N)==F(summary80['H_integral_upper'])
    nodes={};rootkeys=[]
    for root,b in zip(roots,rows):
        k=key(5,root['indices']);rootkeys.append(k)
        nodes[k]={'level':5,'indices':root['indices'],'weight':b['weight'],'selfL':b['weight']*b['integral_lower_i'],'selfU':b['weight']*b['integral_upper_i'],'children':None,'parent':None}
    frontier=set(rootkeys);batch_results=[];alljobs={};allcells={}
    base82=path['082_manager'].parent
    for bid,parents_expected,count in [(0,64,677),(1,256,3102)]:
        tag=f'{bid:04d}';br=read(base82/f'batch{tag}_registration.json');rc=read(base82/f'batch{tag}_receipt.json');jobs=read(base82/f'batch{tag}_jobs.json.gz');cells=read(base82/f'batch{tag}_cells.json.gz')
        assert rc['status']=='complete' and rc['workers_joined'] and len(jobs)==len(cells)==count
        assert sha(base82/f'batch{tag}_jobs.json.gz')==br['jobs_sha256'] and sha(base82/f'batch{tag}_cells.json.gz')==rc['cell_certificates_sha256']
        chosen=sorted(frontier,key=lambda k:(-(nodes[k]['selfU']-nodes[k]['selfL']),nodes[k]['level'],nodes[k]['indices']))[:parents_expected]
        assert [r['key'] for r in br['selected_parents']]==chosen
        assert all(r['gap_units']==nodes[r['key']]['selfU']-nodes[r['key']]['selfL'] for r in br['selected_parents'])
        cm={c['key']:c for c in cells};jm={j['key']:j for j in jobs};assert len(cm)==len(jm)==count and set(cm)==set(jm)
        ordered_count=boundary_count=0
        for pk in chosen:
            parent=nodes[pk];cs=br['children_by_parent'][pk]
            assert set(cs)=={key(parent['level']+1,c) for c in children(tuple(parent['indices']))}
            ordered_count+=partition(tuple(parent['indices']),[idx(k) for k in cs],True);boundary_count+=81
            assert all(k not in nodes for k in cs)
            for k in cs:
                j,b=jm[k],cm[k];assert j['parent']==pk and j['level']==parent['level']+1 and j['indices']==list(idx(k))
                w=weight(tuple(j['indices']));assert w==j['weight']==b['weight']
                geometry(j['point'],j['level'],j['indices']);D,ns=projections(j['point'],vertices)
                mid=[F(2*i+1,2**(j['level']+1)) for i in j['indices']];eta=F(1,2**j['level']*4096)
                candidates=[[F(1)]+mid]+[[F(1)]+[x+eta*beta**n for n,x in enumerate(mid)] for beta in (F(1,8),F(1,16),F(1,32))]
                jits=[(-1,-2,-4,-8),(-1,-3,-7,-11),(-2,-5,-11,-19),(-3,-7,-17,-29),(-2,-11,-23,-43),(-5,-13,-31,-47),(-7,-19,-37,-49)]
                candidates += [[F(1)]+[x+F(v,3351328) for x,v in zip(mid,jit)] for jit in jits]
                assert list(map(F,j['point']))==candidates[j['candidate']] and len(j['qualification'])==j['candidate']+1
                for d,p in zip(j['qualification'],candidates):
                    DD=lcm(*(x.denominator for x in p));nn=[sum(int(x*DD)*z for x,z in zip(p,v)) for v in vertices]
                    inside=all(F(i,2**j['level'])<x<F(i+1,2**j['level']) for i,x in zip(j['indices'],p[1:]))
                    assert d['D']==DD and d['inside']==inside and d['distinct']==len(set(nn)) and d['minimum_absolute_projection']==min(map(abs,nn))
                assert all(not(d['inside'] and d['distinct']==280) for d in j['qualification'][:-1])
                assert j['cached_H1over64'] is None
                chord_check(b['bounds'],j['level'])
                nodes[k]={'level':j['level'],'indices':j['indices'],'weight':w,'selfL':w*b['bounds']['integral_lower_i'],'selfU':w*b['bounds']['integral_upper_i'],'children':None,'parent':pk}
                alljobs[k]=j;allcells[k]=b
            parent['children']=cs;frontier.remove(pk);frontier.update(cs)
        assert frontier==set(read(base82/f'batch{tag}_frontier.json.gz'))
        assert sum((F(nodes[k]['weight'],1<<(4*nodes[k]['level'])) for k in frontier),F(0))==1
        # Different evaluation order: level-descending dynamic programming, not recursive DFS.
        for k in sorted(nodes,key=lambda k:nodes[k]['level'],reverse=True):
            n=nodes[k];cs=n['children'];n['effL']=max(n['selfL'],sum(nodes[c]['effL'] for c in cs)) if cs else n['selfL'];n['effU']=min(n['selfU'],sum(nodes[c]['effU'] for c in cs)) if cs else n['selfU'];assert n['effL']<=n['effU']
        lh=F(sum(nodes[k]['effL'] for k in rootkeys),N);uh=F(sum(nodes[k]['effU'] for k in rootkeys),N)
        a=rc['after'];assert lh==F(a['LH_truncated_integral']) and uh==F(a['UH_truncated_integral'])
        assert F(sum(nodes[k]['selfU'] for k in frontier),N)==F(a['physical_raw_U']) and F(sum(nodes[k]['selfL'] for k in frontier),N)==F(a['physical_raw_L'])
        assert len(frontier)==a['physical_frontier_count'] and F(a['physical_weighted_volume'])==1
        E=F(a['fixed_E076_upper']);B=F(a['B_upper']);I=F(a['I_upper']);g=F(a['conditional_g']);cost=F(a['conditional_cost_upper']);Bl=F(a['fixed_E_budget_lower_diagnostic_B'])
        assert B>E and (B-E)**5>=uh**4 and I**4>=B**5 and cost**4>=(5*B/g)**5
        assert Bl>E and (Bl-E)**5<=lh**4
        assert a['conditional_cost_strictly_below32_exact']==(5**5*B**5<32**4*g**5)
        assert a['fixed_E_budget_route_impossible_below32']==(5**5*Bl**5>=32**4*g**5)
        assert not a['g_independently_qualified_here'] and a['no_tube_added']
        cap=F(ceil(I*10**12),10**12);full_I_cap(uh,E,cap)
        batch_results.append({'batch':bid,'parents':parents_expected,'children':count,'ordered_children_checked':ordered_count,'closed_boundary_points_checked':boundary_count,'frontier_count':len(frontier),'hierarchy_count':len(nodes),'LH':str(lh),'UH':str(uh),'independent_simple_I_cap':str(cap),'fixed_recipe_cost_below32':a['conditional_cost_strictly_below32_exact']})
        print('082 checked batch',bid,'children',count,'seconds',round(time.monotonic()-start,3),flush=True)
    finalhier=read(base82/'numeric_hierarchy.json.gz')
    for k,n in finalhier['expanded_nodes'].items():
        for field in ['level','indices','weight','selfL','selfU','effL','effU','parent','children']:assert nodes[k][field]==n[field]
    assert len(finalhier['expanded_nodes'])==320
    for k,n in finalhier['root_updates'].items():
        for field in ['selfL','selfU','effL','effU','children']:assert nodes[k][field]==n[field]
    abstract=[]
    for m in reg['abstract_m']:
        ps=[tuple(reversed(x)) for x in combinations_with_replacement(range(m),4)]
        assert sum(weight(i) for i in ps)==m**4
        for p in ps:partition(p,children(p))
        abstract.append({'m':m,'canonical_count':len(ps),'weight_sum':m**4})
    synth=synthetic_cuts()
    samples80=read(path['080_samples'])['samples'];samples82=read(HERE/'sample082.json.gz');sample_results=[]
    for source,entries in [('080',samples80),('082',[r['evidence'] for r in samples82])]:
        for ev in entries:
            if source=='080':
                rk=ev['rank'];job=roots[rk];b=rows[rk]
                assert hashlib.sha256(json.dumps(ev,sort_keys=True,separators=(',',':')).encode()).hexdigest()==b['mathematical_sha256']
                level=5;assert ev['point']==job['point'] and ev['endpoints'][24]['H']==job['cached_H1over64']
                assert ev['endpoints'][24]['origin']=='cached066_H1over64'
            else:
                rk=ev['rank'];job=alljobs[rk];b=allcells[rk]['bounds'];level=job['level'];assert ev['point']==job['point'] and ev['bounds']==b
            assert ev['geometry']==geometry(ev['point'],level,job['indices'])
            D,ns=projections(ev['point'],vertices);assert D==ev['D']
            ws={int(k):F(v) for k,v in ev['weights'].items()};assert len(ws)==140 and set(ws)=={k*k for k in ns}
            assert sum(ws.values())==0 and sum(w*k for k,w in ws.items())==0 and sum(w*k*k for k,w in ws.items())==-F(detsum,2)
            ep=ev['endpoints'];assert [e['t'] for e in ep]==read(path['080_registration'])['time_endpoints'] and len(ep)==73
            masses=[]
            for i,e in enumerate(ep):
                lo,hi=F(e['H']['lower']),F(e['H']['upper']);assert 0<lo<=hi
                nxt=ep[i+1]['H'] if i+1<73 else {'lower':'0','upper':'0'}
                ml=lo-F(nxt['upper']);mu=hi-F(nxt['lower']);assert 0<ml<=mu;masses.append((ml,mu))
            assert sum(m[0] for m in masses)<=F(ep[0]['H']['lower'])<=F(ep[0]['H']['upper'])<=sum(m[1] for m in masses)
            sample_results.append({'source':source,'id':rk,'D':D,'positive_bands':73,'weight_moments_passed':True})
    assert len(samples80)==13 and len(samples82)==8
    scalar=read(path['080_cap']);short=read(path['076_results']);E=F(scalar['E076_norm']['upper']);eps=F(1,512);T=F(1,16)
    assert scalar['E076_norm']==next(c['E_norm'] for c in short['cases'] if c['epsilon']=='1/512')
    C3,C1,piL,P5=map(F,[scalar['C3_upper'],scalar['C1_upper'],scalar['pi_lower'],scalar['P5_upper']])
    assert C3==F(short['derived_constants']['C3']['upper']) and C1==F(short['derived_constants']['C1']['upper'])
    assert piL==F(short['derived_constants']['pi_outward_384_dyadic']['lower']) and P5*P5*(1<<33)*piL**6>=1
    assert sum((F(7,10)**j/F(factorial(j)) for j in range(5)),F(0))>2
    roughcap=C3/2*(eps**(-2)-T**(-2))+C1*F(7,2)+P5/(4*T**4);assert roughcap<3
    assert E**5>=F(5,6)**5*eps**6*(C3+C1*eps*eps)/piL**4 and E<F(1,50000)
    full_I_cap(F(cp['upper_units'],N),E,F(summary80['I_upper_simple']))
    ast80=ast.parse(path['080_engine_v2'].read_text());ast82=ast.parse(path['082_kernel'].read_text())
    for fname in ['bounds','qform']:
        a=[ast.dump(n,include_attributes=False) for n in ast80.body if isinstance(n,ast.FunctionDef) and n.name==fname]
        b=[ast.dump(n,include_attributes=False) for n in ast82.body if isinstance(n,ast.FunctionDef) and n.name==fname];assert a==b
    out={'status':'passed','pid':os.getpid(),'source_sha256':sha(__file__),'fingerprint':fp,'elapsed_seconds':time.monotonic()-start,'manifest_file_counts':manifest_counts,'080_complete_roots_and_anchors':52360,'080_positive_projection_qualification_count':52360,'080_weight_sum':cp['weight'],'080_integer_fourth_root_chords':52360,'082_integer_fourth_root_chords':3779,'080_LH':str(F(cp['lower_units'],N)),'080_UH':str(F(cp['upper_units'],N)),'080_simple_I_cap':summary80['I_upper_simple'],'batches':batch_results,'abstract_layers':abstract,'synthetic_33cut_cases':synth,'endpoint_samples':sample_results,'fullface_cap3_rational_upper':str(roughcap),'E076_norm_upper':str(E),'no_new_H_log_exp_Arb_stencil_jet_calls':True,'source080_082_bounds_qform_AST_identical':True,'085_not_consumed':True,'sample_scope':'13 frozen080 samples and8 first082-batch samples, no endpoint reevaluation','new_global_C5_claim':False}
    atomic(HERE/'results.json',out);cp['status']='complete';cp['elapsed_seconds']=out['elapsed_seconds'];atomic(cpfile,cp)
    print(json.dumps({k:out[k] for k in ['status','pid','elapsed_seconds','080_integer_fourth_root_chords','082_integer_fourth_root_chords','batches','no_new_H_log_exp_Arb_stencil_jet_calls']},indent=2))


if __name__=='__main__':main()
