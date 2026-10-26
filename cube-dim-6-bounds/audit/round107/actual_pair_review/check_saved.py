#!/usr/bin/env python3
"""One independent saved/exact107 audit. No producer imports or pair replay."""
import sys
sys.dont_write_bytecode=True
import os,json,hashlib,time
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from datetime import datetime,timezone
OUT=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round107/actual_pair_geometry')
START=time.monotonic();CPU=time.process_time();HASHES={}
def need(x,msg):
    if not x:raise ValueError(msg)
def sha(p):
    p=Path(p);v=hashlib.sha256(p.read_bytes()).hexdigest();old=HASHES.setdefault(str(p),v)
    need(old==v,'Changed evidence '+str(p));return v
def read(p):sha(p);return json.loads(Path(p).read_text())
def save(name,data):
    with (OUT/name).open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
def interval(x):return Q(x['lower']),Q(x['upper'])
def dyadic(x):
    s=Q(2)**int(x['common_binary_exponent']);return int(x['lower_mantissa'])*s,int(x['upper_mantissa'])*s
def qr(p,a,b):
    pts={a,b}
    if p[2] and a < -p[1]/(2*p[2]) < b:pts.add(-p[1]/(2*p[2]))
    vals={x:p[0]+p[1]*x+p[2]*x*x for x in pts};lo=min(vals.values());hi=max(vals.values())
    return lo,hi,{x for x,v in vals.items() if v==lo},{x for x,v in vals.items() if v==hi}
def rectangle_min(f,ua,ub,aa,ab):
    # Exact quadratic extremum theorem: corners, edge stationary points,
    # interior gradient zero. A singular stationary line meets a boundary.
    c=f(Q(0),Q(0));u=(f(Q(1),Q(0))-f(Q(-1),Q(0)))/2
    a=(f(Q(0),Q(1))-f(Q(0),Q(-1)))/2
    uu=(f(Q(1),Q(0))+f(Q(-1),Q(0)))/2-c
    av=(f(Q(0),Q(1))+f(Q(0),Q(-1)))/2-c
    cross=f(Q(1),Q(1))-c-u-a-uu-av
    pts={(x,y) for x in (ua,ub) for y in (aa,ab)}
    for x in (ua,ub):
        if av:
            y=-(a+cross*x)/(2*av)
            if aa<y<ab:pts.add((x,y))
    for y in (aa,ab):
        if uu:
            x=-(u+cross*y)/(2*uu)
            if ua<x<ub:pts.add((x,y))
    det=4*uu*av-cross*cross
    if det:
        x=(-2*av*u+cross*a)/det;y=(-2*uu*a+cross*u)/det
        if ua<x<ub and aa<y<ab:pts.add((x,y))
    return min(f(x,y) for x,y in pts)
def main():
    reg=read(OUT/'registration.json')
    need(sha(Path(__file__))==reg['checker_sha256'],'Registered checker')
    need(sha(AUTHOR/'manifest.json')==reg['author_manifest_sha256'],'Frozen author seal')
    manifest=read(AUTHOR/'manifest.json')
    for row in manifest['files']:
        p=AUTHOR/row['name'];need(p.resolve().is_relative_to(AUTHOR.resolve()),'Author whitelist path')
        need(sha(p)==row['sha256'] and p.stat().st_size==row['bytes'],'Author saved bytes')
    for row in read(AUTHOR/'source_manifest.json')['sources']:need(sha(row['path'])==row['sha256'],'Frozen scientific/admin source')
    administrative=read(AUTHOR/'administrative_source_snapshots.json')
    need(sha(AUTHOR/'administrative_source_snapshots.json')==reg['administrative_asof_sha256'] and administrative['files']==[],'107 uses immutable amendment directly; no mutable source substitution')
    data=read(AUTHOR/'results.json');geo=read(AUTHOR/'pair_geometry_certificate.json');cells=read(AUTHOR/'cell_certificate.json')
    up=read(AUTHOR/'upstream_constants.json');summary=read(AUTHOR/'exact_summary.json')
    release=read(AUTHOR/'parent_release.json');receipt=read(AUTHOR/'run_receipt.json');tool=read(AUTHOR/'tool_observation.json')
    need(release['authorized'] and release['script_sha256']==sha(AUTHOR/'run_prepared.py')==data['script_sha256'] and
         release['registration_sha256']==sha(AUTHOR/'registration.json')==data['registration_sha256'],'Original released source')
    need(receipt['status']=='COMPLETED' and receipt['pid']==52476 and receipt['completed_configurations']==2 and receipt['failure'] is None and
         receipt['retries']==receipt['children']==0 and receipt['release_sha256']==sha(AUTHOR/'parent_release.json'),'Original actual receipt')
    need(tool['scientific_invocations']==1 and tool['initial_exec_command']['session_id']==1726 and
         tool['terminal_write_stdin']['exit_code']==0 and tool['terminal_write_stdin']['pid']==52476,'Original handle/exit binding')
    cfg={x['configuration']:x for x in data['configurations']};c=cfg[2];delta=Q(1,64)
    need(data['Cn_improved'] is False and data['105_or108_unaccepted_bound_applied_as_theorem'] is False,'107 labels remain conditional')
    need(geo['k']==cells['k']==1 and cells['delta']=='1/64' and cells['q']=='6/5','Exact107 domain')
    # Independent tail-intersection combinatorics, never iterate1860 body pairs.
    tables={};checks=0
    for branch,excluded in [('v',5),('w',0)]:
        table={}
        for a in range(6):
            for b in range(6):
                if a==excluded or b==excluded:continue
                for t in range(max(0,a+b-5),min(a,b)+1):
                    h=a+b-2*t
                    if not h:continue
                    n=comb(5,a)*comb(a,t)*comb(5-a,b-t);dot=6-2*h;D=36-dot*dot
                    A=2*a-5;B=2*b-5
                    p=(Q(1)-Q(12-2*dot,D),-Q((12-2*dot)*(A+B),D),Q(5)-Q(6*A*A-2*dot*A*B+6*B*B,D))
                    table[p,D]=table.get((p,D),0)+n;checks+=1
        saved=geo['branches'][branch];rows=saved['classes'];actual={(tuple(map(Q,x['P_coefficients'])),x['J2_squared']):x['ordered_multiplicity'] for x in rows}
        need(table==actual and len(rows)==len(actual)==23 and sum(table.values())==saved['pair_visits']==930,'Independent23classes/full930 '+branch)
        need({D:sum(n for (p,d),n in table.items() if d==D) for D in (20,32,36)}=={20:180,32:450,36:300},'Determinants/multiplicities')
        for row in rows:
            lo,hi,lpts,hpts=qr(tuple(map(Q,row['P_coefficients'])),Q(1,4),Q(3,4));gl,gu=dyadic(row['normal_radius_gap_at_full_delta'])
            need(lo==Q(row['P_min_exact']) and hi==Q(row['P_max_exact']) and lpts==set(map(Q,row['minimum_nodes'])) and hpts==set(map(Q,row['maximum_nodes'])),'Full-domain exact quadratic extrema')
            need(lo>delta**2 and 0<gl<=gu and (gl+delta)**2<=lo<=(gu+delta)**2,'Exact rational individual root-gap enclosure')
        tables[branch]=table
    # The first-line inequalities are verified as quadratic nonnegativity on
    # their entire two-parameter rectangles, not by testing source points.
    for branch,cap,excluded in [('v',Q(1,2),5),('w',Q(1,16),-5)]:
        for S in (-5,-3,-1,1,3,5):
            if S==excluded:continue
            pieces=[(Q(1,4),Q(1,2),lambda u:u),(Q(1,2),Q(3,4),lambda u:1-u)] if branch=='v' else [(Q(1,4),Q(3,4),lambda u:Q(7,8)-u)]
            for left,right,m in pieces:
                def f(u,a):
                    A=1-a;B=u-a if branch=='v' else u+a
                    return 6*A*A+30*B*B-(A+S*B)**2-5*m(u)**2
                need(rectangle_min(f,left,right,-cap,cap)>=0,'All-source first-line gap lemma')
    from flint import arb,fmpq,ctx
    ctx.prec=512
    def A(x):
        x=Q(x);return arb(fmpq(x.numerator,x.denominator))
    def ep(x):
        n,e=x.lower().man_exp();m,f=x.upper().man_exp();return Q(int(n))*Q(2)**int(e),Q(int(m))*Q(2)**int(f)
    def inside(x,saved,msg):
        l,r=interval(saved);a,b=ep(x);need(l<=a<=b<=r,msg)
    pi=arb.pi();ap=1/(2*A(5).sqrt());mass=A(Q(10,7))*pi*pi/2*A(5).sqrt()*A(delta)**A(Q(14,5))
    inside(ap,cells['large_shape_coefficient_1_over_Jvw'],'Large branch coarea')
    inside(mass,cells['mass_constant_after_s13_exact_integration'],'s13/Jac mass coefficient')
    masslo,masshi=interval(cells['mass_constant_after_s13_exact_integration'])
    # Recover the tree from seeds and every recorded split; a volume sum alone
    # is insufficient. Every leaf must match its full reconstructed rectangle.
    ub=list(map(Q,cells['u_clamp_partition']));need(ub==[Q(1,4),Q(11,40),Q(2,5),Q(1,2),Q(3,4)],'Clamp/line partition')
    expected={(a,b,s,t) for a,b in zip(ub,ub[1:]) for s,t in [(Q(0),Q(1,2)),(Q(1,2),Q(1))]}
    live={};seen=set()
    for row in cells['initial_seeds']:
        rect=tuple(map(Q,row['u']+row['s']));cid=row['id'];need(rect in expected and cid not in seen,'Exact seed tiling')
        expected.remove(rect);live[cid]=(*rect,0,0,cid);seen.add(cid)
    need(not expected and len(live)==8,'All8 seeds')
    for split in cells['refinements']:
        p=split['parent'];need(p in live,'Refine an existing leaf exactly once');a,b,s,t,du,ds,seed=live.pop(p);ids=split['children']
        need(len(ids)==2 and len(set(ids))==2 and all(x not in seen for x in ids),'Fresh two-child tree')
        if split['axis']=='u':
            need(du<=ds,'Registered alternating axis');mid=(a+b)/2;rects=[(a,mid,s,t,du+1,ds,seed),(mid,b,s,t,du+1,ds,seed)]
        else:
            need(split['axis']=='s' and du>ds,'Registered alternating axis');mid=(s+t)/2;rects=[(a,b,s,mid,du,ds+1,seed),(a,b,mid,t,du,ds+1,seed)]
        for cid,rect in zip(ids,rects):live[cid]=rect;seen.add(cid)
    leaves=cells['leaves'];need(len(leaves)==len(live)==1024 and len(seen)==2040 and len(cells['refinements'])==1016,'1024leaf/2040node complete tree')
    need({x['id'] for x in leaves}==set(live),'Exact final-leaf set')
    rootD={D:A(D).sqrt() for D in (20,32,36)};root56=A(Q(5,6)).sqrt();cache={};sumlo=Q(0);sumhi=Q(0);mingap=None
    for row in leaves:
        a,b,s,t,du,ds,seed=live[row['id']];rlo=delta*s**5;rhi=delta*t**5
        need(tuple(map(Q,row['u']+row['s']))==(a,b,s,t) and row['depths']==[du,ds] and row['seed']==seed and du+ds<=20,'Leaf rectangle/depth/seed')
        need(list(map(Q,row['r']))==[rlo,rhi],'Exact desingularized radius endpoints')
        g=Q(row['positive_denominator_lower']);need(g>0,'Positive recorded cell denominators');mingap=g if mingap is None else min(mingap,g)
        low=A(0);high=A(0)
        for branch in ('v','w'):
            Wlo=A(0);Whi=A(0)
            for (p,D),n in tables[branch].items():
                key=(p,a,b)
                if key not in cache:
                    pl,pu,_,_=qr(p,a,b);cache[key]=(pl,pu,A(pl).sqrt(),A(pu).sqrt())
                pl,pu,rl,ru=cache[key]
                need(pl>=(g+rhi)**2,'Every saved pair denominator >=recordedg')
                Wlo+=A(n)/(rootD[D]*(ru-A(rlo)));Whi+=A(n)/(rootD[D]*(rl-A(rhi)))
            if branch=='v':
                p=(Q(5,6),-Q(5,3),Q(5,6)) if b<=Q(2,5) else (Q(3,2),-Q(5),Q(5))
                need(b<=Q(2,5) or a>=Q(2,5),'No v clamp crossing')
                mlo,mhi=(a,b) if b<=Q(1,2) else (1-b,1-a)
                need(b<=Q(1,2) or a>=Q(1,2),'No line switch crossing')
            else:
                p=(Q(5,6),Q(5,3),Q(5,6)) if b<=Q(11,40) else (Q(147,128),-Q(5,8),Q(5))
                need(b<=Q(11,40) or a>=Q(11,40),'No w clamp crossing');mlo,mhi=Q(7,8)-b,Q(7,8)-a
            flo,fhi,_,_=qr(p,a,b)
            need(flo>=(g+rhi)**2 and Q(5,6)*mlo*mlo>=(g+rhi)**2,'Every saved origin/line denominator >=recordedg')
            h0lo=A(flo).sqrt()-A(rhi);h0hi=A(fhi).sqrt()-A(rlo)
            h1lo=root56*A(mlo)-A(rhi);h1hi=root56*A(mhi)-A(rlo)
            low+=Wlo/(h0hi*h0hi*h1hi*h1hi);high+=Whi/(h0lo*h0lo*h1lo*h1lo)
        tl,tu=interval(row['regular_shape_T']);need(0<tl<=ep(ap+A(rlo)*low)[0] and ep(ap+A(rhi)*high)[1]<=tu,'Full-cell saved regular T covers independent512 geometry')
        pl,pu=interval(row['regular_shape_q_power']);need(0<pl<=pu and pl**5<=tl**6 and pu**5>=tu**6,'Exact fractional-power rational inequalities')
        weight=(b-a)*(t**14-s**14);need(weight==Q(row['exact_weight'])>0,'Integrated s13 weight')
        lo,hi=interval(row['weighted_shape_integral']);need(lo==masslo*weight*pl and hi==masshi*weight*pu,'Complete rational leaf endpoints')
        sumlo+=lo;sumhi+=hi
    need((sumlo,sumhi)==interval(cells['shape_integral'])==interval(c['shape_integral']),'Exact1024leaf rational endpoint sum')
    need(c['leaves']==1024 and c['cell_bounds']==2040 and c['stopping_outcome']=='FINITE_NODE_CAP_REACHED','Actual fixed cap')
    Kold=3*A(3).sqrt()/(64*pi**4);Kref=1/(8*pi**5);q=A(Q(6,5))
    for value,key in [(Kold,'Kold'),(Kref,'Kref_conditional_only'),(Kold**q,'scale_old'),(Kref**q,'scale_ref_conditional_only')]:
        l,r=dyadic(c[key]);a,b=ep(value);need(l<=a<=b<=r,'Original coefficient/scale interval '+key)
    for label,key in [('old098_valid','scale_old'),('ref108_conditional_shape_only','scale_ref_conditional_only')]:
        l,r=dyadic(c[key]);need(interval(c['benchmarks'][label]['I'])==(l*sumlo,r*sumhi),'Same one-shape coefficient scaling')
    need(c['benchmarks']['ref108_conditional_shape_only']['coefficient_validity']=='108_NOT_ACCEPTED_CONDITIONAL_SHAPE_EXPRESSION_ONLY' and c['108_theorem_used'] is False,'Do not relabel pending108')
    up104=read(AUTHOR.parents[1]/'round104/complement_event_budget/results.json');need(up['old104_k1_weighted_I']==up104['configurations'][2]['cases'][0]['weighted_I'],'Identical saved104 k1 object')
    paid062=read('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round062/primary_joint_account/actual_run/results.json')
    need(Q(up['paid_CEH_I_upper'])==Q(paid062['paid_C_E_H_I_upper']),'Frozen actual fullG paidCEH source')
    oldlo,oldhi=interval(c['benchmarks']['old098_valid']['I']);need(oldhi<=Q(2262318643212952363,10**24)<dyadic(up['old104_k1_weighted_I'])[0],'Strict same-region Kold upper improvement')
    target=64*A(Q(1,390000))**q;paid=Q(up['paid_CEH_I_upper']);rl,ru=interval(c['conditional_remaining_I_interval']);ta,tb=ep(target)
    need(rl<=ta-paid<=tb-paid<=ru,'Conditional remaining target interval')
    for label in ('old098_valid','ref108_conditional_shape_only'):
        lo,hi=interval(c['benchmarks'][label]['I']);need(interval(c['benchmarks'][label]['S_plan_other_regions_zero'])==(rl-hi,ru-lo),'Conditional remaining endpoint subtraction')
    mapping={'shape_expression':interval(c['shape_integral']),'old098_valid_expression':(oldlo,oldhi),
             'old098_partial_remaining_allowance':interval(c['benchmarks']['old098_valid']['S_plan_other_regions_zero']),
             'ref108_unaccepted_conditional_expression':interval(c['benchmarks']['ref108_conditional_shape_only']['I']),
             'ref108_unaccepted_partial_remaining_allowance':interval(c['benchmarks']['ref108_conditional_shape_only']['S_plan_other_regions_zero'])}
    for label,(lo,hi) in mapping.items():
        item=summary[label];need(Q(item['lower'])==Q(item['lower_decimal'])<=lo<=hi<=Q(item['upper_decimal'])==Q(item['upper']),'Complete outward rational decimal summary')
    need(summary['results_sha256']==sha(AUTHOR/'results.json') and summary['actual_I_lower_supplied'] is False and summary['complete_global_S'] is False,'Partial-I scope')
    for p,v in list(HASHES.items()):need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v,'Postcheck original mutation')
    result=dict(status='PASS_FULL107_SAVED_CERTIFICATE_INDEPENDENT_GEOMETRY_AND_ALL_LEAF_ENDPOINTS',source_hashes=HASHES,
        author_manifest_sha256=sha(AUTHOR/'manifest.json'),hypergeometric_cases=checks,classes_each_branch=23,ordered_pairs_each_branch=930,
        actual_pair_iterations=0,geometry_range_precision_bits=512,all_leaf_endpoints_checked=1024,complete_tree_nodes=2040,
        exact_shape_interval=list(map(str,(sumlo,sumhi))),Kold_k1_upper=str(Q(2262318643212952363,10**24)),
        minimum_recorded_cell_gap=str(mingap),clamp_and_all_source_line_gap_exact=True,s13_and_sqrt5_Jacobian_review=True,
        pending108_label_unchanged=True,global_g_otherk_far_fullI_S_unpaid=True,true_I_lower_or_C6_lower=False,Cn_improved=False,
        producer_replays=0,old_kernels=0,actual_native_kernel_calls=0,MATLAB=0,source_writes=0,children_created=0)
    save('checked_saved.json',result)
    save('run_receipt.json',dict(status='PASS',pid=os.getpid(),utc=datetime.now(timezone.utc).isoformat(),
        wall_before_receipt_after_result_save=time.monotonic()-START,cpu_before_receipt=time.process_time()-CPU,
        runtime_limit=None,first_failure_stop=True,attempts=1,author_generator_runs=0,old_kernel_native_calls=0,
        result_sha256=sha(OUT/'checked_saved.json')))
    print(json.dumps(dict(status=result['status'],leaves=1024,classes_each_branch=23,Kold_k1_upper=result['Kold_k1_upper'],Cn_improved=False)))
if __name__=='__main__':
    try:main()
    except BaseException as error:
        save('failure_receipt.json',dict(status='STOP_NO_RETRY',pid=os.getpid(),error=repr(error),wall=time.monotonic()-START,cpu=time.process_time()-CPU,producer_replays=0))
        raise
