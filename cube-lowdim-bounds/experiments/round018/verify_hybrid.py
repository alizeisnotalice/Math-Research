#!/usr/bin/env python3
"""Independent replay interface for a saved hybrid attempt; no adaptive work.

Default compact-summary mode checks auxiliary/cap evidence and at most three
saved numeric leaves. --full FILE additionally checks the full partition tree;
if that attempt is closed, every numeric leaf is recomputed. Reuse only the
reviewed round008 Fraction/log kernel and round017 Poly/AD primitives. Never
import hybrid_certify.py or call its evaluate/complete_bundle/main functions.
"""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import argparse
import importlib.util
import json
import sys
from time import monotonic

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
R=Q(15,8)
E={2:Q(69,100),6:Q(19,50),7:Q(11,10),8:Q(3,4)}
K=(Q(99,100),Q(1),Q(99,100),Q(1))
ROOT_BOX=(Q(1,100),Q(1),Q(0),Q(1))


def need(value,label):
    if not value:
        raise ValueError(label)


def box(row):
    return tuple(Q(row[k]) for k in ('r_lo','r_hi','q_lo','q_hi'))


def digest(data):
    return sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=obj
    spec.loader.exec_module(obj)
    return obj


def resolve(name,wanted,mapping):
    candidates=[]
    if name in mapping:
        p=Path(mapping[name]);candidates.append(p if p.is_absolute() else ROOT/p)
    if name=='round017/boundary_cap_verification.json':
        candidates.append(ROOT/'cap_input_record.json')
    candidates.append(ROOT.parent/name)
    for p in candidates:
        if p.is_file() and sha256(p.read_bytes()).hexdigest()==wanted:
            return p
    raise FileNotFoundError('No matching dependency bytes: '+name+' (use --dependency-map).')


def prepare_primitives(data,mapping):
    core_path=resolve(data['core_dependency'],data['core_sha256'],mapping)
    base_name='round017/compare_bounds.py'
    base_path=resolve(base_name,data['dependencies'][base_name],mapping)
    c=load('round018_replay_fraction_log',core_path)
    base=load('round018_replay_reviewed_primitives',base_path)
    c.R=R
    base.AD.I,base.AD.log_i=c.I,staticmethod(c.log_i)
    # Rebuild the formula here; do not call the certification/old bound builders.
    r,q=base.Poly({(1,0):1}),base.Poly({(0,1):1})
    u,v,s=r*(1+q)/2,r*(1-q)/2,1-r
    a,b,zp,zm=1-v,1-u,1+r*q,1-r*q
    polynomials={}
    for p in E:
        n=p-1
        P=sum((a**(n-k)*(zp**k+s**k)+b**(n-k)*(zm**k+s**k))/Q(n-k) for k in range(n))
        polynomials[p]=(P,zp**n+s**n,zm**n+s**n)
    T=2/(R-1)*(2*sum(r**j/Q((j+2)*R**(j+1)) for j in range(29))+2*sum((r*q)**(2*k)/Q((2*k+2)*R**(2*k+1)) for k in range(15)))
    derivatives=[T,T.diff(0),T.diff(1),T.diff(0).diff(0),T.diff(0).diff(1),T.diff(1).diff(1)]
    need(all(all(x>=0 for x in P.d.values()) for P in derivatives),'positive tail derivative coefficients')
    return c,base,polynomials,derivatives


class NoAdaptiveGuard:
    def check(self):
        pass


def bundles(c,base,B,polynomials,tail,full):
    rl,rh,ql,qh=B
    rc,qc=(rl+rh)/2,(ql+qh)/2
    AD=base.AD
    rp,qp=AD.variable(c.I.point(rc),0),AD.variable(c.I.point(qc),1)
    rb,qb=AD.variable(c.I(rl,rh),0),AD.variable(c.I(ql,qh),1)
    def variables(r,q):
        return r*(1+q)/2,r*(1-q)/2,1-r,1+r*q,1-r*q
    up,vp,sp,zpp,zmp=variables(rp,qp)
    ub,vb,sb,zpb,zmb=variables(rb,qb)
    if full:
        need(all(x.v.lo>0 for x in (up,vp,sp,zpp,zmp,ub,vb,sb,zpb,zmb)),'full replay requires strictly positive arguments')
    guard=NoAdaptiveGuard()
    point,uniform=base.tail_bundles(tail,B,guard)
    for p,(P,Cu,Cv) in polynomials.items():
        Pp,Pb=base.shifted_bundles(P,B,guard)
        Cp,Cb=base.shifted_bundles(Cu,B,guard)
        n,alpha=p-1,2*p*E[p]
        if full:
            Dp,Db=base.shifted_bundles(Cv,B,guard)
            point-=alpha*(Pp+Cp*up.log()+Dp*vp.log()-zpp**n*zpp.log()-zmp**n*zmp.log()-2*sp**n*sp.log())
            uniform-=alpha*(Pb+Cb*ub.log()+Db*vb.log()-zpb**n*zpb.log()-zmb**n*zmb.log()-2*sb**n*sb.log())
        else:
            point-=alpha*(Pp+Cp*up.log()-zpp**n*zpp.log())
            uniform-=alpha*(Pb+Cb*ub.log()-zpb**n*zpb.log())
    return point,uniform


def taylor(base,point,uniform,B):
    rr,rq=(B[1]-B[0])/2,(B[3]-B[2])/2
    mean=point.v.lo-base.maxabs(uniform.g[0])*rr-base.maxabs(uniform.g[1])*rq
    second=point.v.lo-base.maxabs(point.g[0])*rr-base.maxabs(point.g[1])*rq-(base.maxabs(uniform.h[0][0])*rr**2+2*base.maxabs(uniform.h[0][1])*rr*rq+base.maxabs(uniform.h[1][1])*rq**2)/2
    record={'center_value':point.v.pair(),'center_gradient':[x.pair() for x in point.g],
            'uniform_gradient':[x.pair() for x in uniform.g],
            'uniform_hessian':[[x.pair() for x in row] for row in uniform.h],
            'mean_lower':str(mean),'second_order_lower':str(second),'best_lower':str(max(mean,second))}
    return max(mean,second),record


def singular(c,B):
    rl,rh,ql,qh=B
    sl,sh,zl,zh=1-rh,1-rl,1-rh*qh,1-rl*ql
    vh=rh*(1-ql)/2
    def w(l,h,n):
        if h==0:return Q(0)
        if l==0:return h**n*c.log_point(h).lo if h<=Q(1,4) else -Q(1,2*n)
        return min(x**n*y for x in (l,h) for y in (c.log_point(l).lo,c.log_point(h).hi))
    C=Q(0);ws={};wz={}
    for p,e in E.items():
        n,alpha=p-1,2*p*e
        C+=alpha*(zl**n+sl**n)
        ws[str(p)],wz[str(p)]=w(sl,sh,n),w(zl,zh,n)
    favorable=-C*min(Q(0),c.log_point(vh).hi) if C and vh else Q(0)
    value=favorable+sum((2*p*e*(wz[str(p)]+2*ws[str(p)]) for p,e in E.items()),Q(0))
    return value,{'C_lower':str(C),'v_upper':str(vh),'favorable_negative_log_lower':str(favorable),
                  'W_s_lower_by_power':{p:str(x) for p,x in ws.items()},
                  'W_zminus_lower_by_power':{p:str(x) for p,x in wz.items()},'total_lower':str(value)}


def numeric_leaf(c,base,polynomials,tail,leaf,allow_minimal=False):
    B=box(leaf)
    need(Q(1,100)<=B[0]<B[1]<=1 and 0<=B[2]<B[3]<=1,'numeric leaf domain')
    evidence=leaf.get('evidence')
    need(evidence is not None or allow_minimal,'detailed evidence or declared minimal format')
    gate={'s':1-B[1],'v':B[0]*(1-B[3])/2,'zminus':1-B[1]*B[3]}
    need(all(x>=0 for x in gate.values()),'nonnegative regularity gate')
    if evidence is not None:
        need({k:str(x) for k,x in gate.items()}==evidence['regularity_lower_endpoints'],'saved gate')
        need(evidence['center']==[str((B[0]+B[1])/2),str((B[2]+B[3])/2)] and evidence['halfwidth']==[str((B[1]-B[0])/2),str((B[3]-B[2])/2)],'center and halfwidth')
    strict=all(x>0 for x in gate.values())
    full_lower=None
    if strict:
        p,u=bundles(c,base,B,polynomials,tail,True)
        full_lower,record=taylor(base,p,u,B)
        if evidence is not None:
            need(record==evidence['full_smooth'],'independent full smooth bundle exactly reproduced')
        if full_lower>=0:
            need(leaf['method']=='full_smooth' and (evidence is None or 'split_boundary' not in evidence),'full acceptance before fallback')
            value=full_lower
        else:
            value=None
    else:
        need(evidence is None or 'full_smooth' not in evidence,'boundary gate forbids full evaluation')
        value=None
    if value is None:
        p,u=bundles(c,base,B,polynomials,tail,False)
        core_lower,record=taylor(base,p,u,B)
        extra,extra_record=singular(c,B)
        split=core_lower+extra
        record.update({'singular_remainder':extra_record,'total_lower':str(split)})
        if evidence is not None:
            need(record==evidence['split_boundary'],'independent split bundle exactly reproduced')
        value=split if full_lower is None else max(split,full_lower)
        expected=('split_boundary' if not strict else 'split_fallback') if full_lower is None or split>=full_lower else 'full_smooth_after_split_fallback'
        need(leaf['method']==expected,'saved route/fallback selection')
    need(value>=0 and str(value)==leaf['generator_lower'],'independently reproduced nonnegative accepted lower')
    return {'node_id':leaf['node_id'],'method':leaf['method'],'generator_lower_display':float(value),'all_saved_bundle_fields_exactly_reproduced':evidence is not None,'exact_saved_lower_reproduced':True}


def auxiliary_and_cap(c,data,mapping):
    aux=data['auxiliary'];cap=data['analytic_cap']
    harmonic=sum((Q(1,k) for k in range(1,8)),Q(0))
    B=2*Q(101,100)**6*harmonic+c.log_point(Q(1,101))+Q(99,100)**7*c.log_point(Q(1,99))
    need(B.pair()==aux['small_radius']['B_interval'] and B.hi<0,'small radius exact interval')
    l2=c.log_point(2)
    H={2:8-8*l2,6:Q(1334,5)-384*l2,7:Q(3108,5)-896*l2,8:Q(149104,105)-2048*l2}
    D=2*(2*c.tail_series(c.I.point(1))+c.I.point(1)/R)/(R-1)
    exterior=D-sum((E[p]*H[p] for p in E),c.I.point(0))
    need({str(p):x.pair() for p,x in H.items()}==aux['exterior']['Hp_intervals'] and D.pair()==aux['exterior']['D_R1_interval'],'exterior constants')
    need(exterior.pair()==aux['exterior']['new_slack_interval'] and exterior.lo>0,'strict exterior condition')
    dominance=R-(R-1)*sum(E.values())/4
    cost=R*R/(R-1)-sum((2*e/Q((p+1)*(p+2)) for p,e in E.items()),Q(0))
    need(dominance==Q(aux['kernel_dominance'])==Q(989,800),'kernel dominance')
    need(cost==Q(aux['cost'])==Q(4841,1260) and Q(aux['previous_strict_upper'])-cost==Q(aux['improvement'])==Q(28363,1118250),'exact candidate mass/cost')
    proof=resolve(cap['proof_dependency'],cap['proof_sha256'],mapping)
    verifier=resolve('round017/verify_boundary_cap.py',cap['verifier_sha256'],mapping)
    saved_path=resolve('round017/boundary_cap_verification.json',cap['saved_verification_sha256'],mapping)
    saved=json.loads(saved_path.read_text())
    need(saved['proof_sha256']==sha256(proof.read_bytes()).hexdigest() and saved['verifier_sha256']==sha256(verifier.read_bytes()).hexdigest(),'cap source evidence chain')
    need(saved['core_sha256']==cap['core_sha256']==data['core_sha256'],'cap reviewed core')
    eps,tau,bmax=Q(2,199),Q(197,199),Q(299,20000)
    omega={2:4*eps*(c.log_point(199).hi-1)};M={2:Q(4)}
    for p in (6,7,8):
        m,n=p-2,p-1
        C=1/(tau*(p-1))+sum((eps**k/Q(m-k) for k in range(m)),Q(0))+eps**m*c.log_point(1/eps).hi
        omega[p]=eps*max(Q(0),2*p*(p-1)*C-4*p)
        M[p]=2*p*(sum((Q(2**k+1,n-k) for k in range(n)),Q(0))+Q(1,2*n))
    U={p:H[p].hi+omega[p]+bmax**(p-1)*M[p] for p in E}
    F=sum((Q(99,100)**j/((j+2)*R**(j+1)) for j in range(29)),Q(0))
    T=2*(2*F+1/R)/(R-1)
    margin=T-sum((E[p]*U[p] for p in E),Q(0))
    need(cap['target']=='real_G' and box(cap['domain'])==K,'cap target and exact domain')
    need(str(T)==cap['tail_lower'] and {str(p):str(x) for p,x in U.items()}==cap['bump_upper_by_power'],'cap tail and bump bounds')
    need(margin>0 and margin==Q(cap['generator_lower'])==Q(saved['strict_cap_generator_margin_fraction']),'cap exact positive margin')
    return {'auxiliary_exact_intervals_reproduced':True,'cap_exact_margin_reproduced':True,
            'cap_margin_display':float(margin),'cap_input_sha256':sha256(saved_path.read_bytes()).hexdigest(),
            'candidate_cost_only':str(cost),'rigorous_upper_before_closed_replay':'411881/106500'}


def tree_check(data,allow_minimal=False):
    nodes,leaves,queue=data['partition_tree'],data['accepted_leaf_data'],data['unresolved_queue']
    need(box(nodes[0])==ROOT_BOX and nodes[0]['parent'] is None,'root domain')
    seen=set();todo=[0];terminals=[];seed=mid=0
    while todo:
        i=todo.pop();node=nodes[i]
        need(i not in seen and node['id']==i,'node reachability/identity')
        seen.add(i);l,h,a,b=box(node)
        need(Q(1,100)<=l<h<=1 and 0<=a<b<=1,'positive tree domain')
        if node['status']!='split':terminals.append(i);continue
        cut=Q(node['cutpoint']);axis=node['axis']
        if node['method']=='seed_partition':
            seed+=1
            need((i,axis,cut) in [(0,'r',Q(99,100)),(2,'q',Q(99,100))],'two exact fixed seed cuts')
        else:
            mid+=1
            need(node['method']=='exactmidpoint','ordinary split method')
            need(cut==((l+h)/2 if axis=='r' else (a+b)/2),'ordinary exact midpoint')
            need(axis==('r' if (h-l)/Q(99,100)>=b-a else 'q'),'normalized longest side')
        need(axis in ('r','q'),'split axis')
        wanted=[(l,cut,a,b),(cut,h,a,b)] if axis=='r' else [(l,h,a,cut),(l,h,cut,b)]
        ids=node['children'];need(len(ids)==2 and ids[0]!=ids[1],'two distinct children')
        for child,B in zip(ids,wanted):
            need(nodes[child]['parent']==i and box(nodes[child])==B,'exact binary parent/child coverage')
            todo.append(child)
    need(len(seen)==len(nodes)==1+2*(seed+mid) and seed==2,'all nodes reachable and seed/midpoint counts')
    accepted={x['node_id'] for x in leaves};pending={x['node_id'] for x in queue}
    need(len(accepted)==len(leaves) and len(pending)==len(queue) and not accepted&pending and accepted|pending==set(terminals),'terminal lists form disjoint exhaustive cover')
    for leaf in leaves:
        node=nodes[leaf['node_id']]
        need(node['status']=='accepted' and node['method']==leaf['method'] and box(node)==box(leaf) and Q(leaf['generator_lower'])>=0,'saved accepted identity/sign')
    for leaf in queue:need(box(nodes[leaf['node_id']])==box(leaf),'unresolved identity')
    caps=[x for x in leaves if x['method']=='analytic_cap']
    need(len(caps)==1 and caps[0]['node_id']==4 and box(caps[0])==K and caps[0]['generator_lower']==data['analytic_cap']['generator_lower'],'one whole audited cap leaf')
    need(caps[0].get('cap_evidence')==data['analytic_cap'] or (allow_minimal and 'cap_evidence' not in caps[0]),'cap detailed evidence or declared minimal format')
    counts=data['counts'];numeric=[x for x in leaves if x['method']!='analytic_cap']
    floors=sum(nodes[i]['status']=='unresolved_width_floor' for i in pending)
    interrupted=sum(nodes[i]['status']=='pending_interrupted_evaluation' for i in pending)
    need(counts['seed_partition_splits']==2 and counts['midpoint_splits']==mid and counts['numeric_accepted_leaves']==len(numeric),'saved split/accept counts')
    need(counts['boxes_completed']==len(numeric)+mid+floors and counts['boxes_interrupted']==interrupted and counts['boxes_started']==counts['boxes_completed']+interrupted<=30000,'started/completed/interrupted accounting')
    need(counts['center_bundle_calls']==counts['uniform_derivative_bundle_calls']==counts['full_smooth_bundle_calls']+counts['split_bundle_calls'],'started bundle accounting')
    need(counts['split_fallback_calls']<=counts['split_bundle_calls']<=counts['boxes_started'] and counts['full_smooth_bundle_calls']<=counts['boxes_started'],'method call bounds')
    need(counts['analytic_cap_leaves']==counts['cap_constant_checks']==1 and counts['LPs']==counts['extra_probe_boxes']==0,'cap and no-extra-work counts')
    need(digest(leaves)==data['accepted_leaf_sha256'] and digest(nodes)==data['partition_tree_sha256'],'full leaf/tree digests')
    def area(rows):return sum(((B[1]-B[0])*(B[3]-B[2]) for B in map(box,rows)),Q(0))
    accepted_area,pending_area=area(leaves),area(queue)
    need(accepted_area+pending_area==Q(99,100),'exact area ledger')
    ledger=data['parameter_area_ledger_only']
    need(Q(ledger['accepted_area'])==accepted_area and Q(ledger['unresolved_area'])==pending_area and Q(ledger['root_area'])==Q(99,100),'saved area ledger')
    closed=not queue and all(x['status'] in ('accepted','split') for x in nodes)
    need(closed==data['complete_partition_cover'],'closed cover flag')
    if closed:
        evaluated_nodes=[n for n in nodes if n.get('method')=='exactmidpoint' or (n['status']=='accepted' and n.get('method')!='analytic_cap') or n['status']=='unresolved_width_floor']
        strict_count=sum(1-box(n)[1]>0 and box(n)[0]*(1-box(n)[3])/2>0 and 1-box(n)[1]*box(n)[3]>0 for n in evaluated_nodes)
        boundary_count=len(evaluated_nodes)-strict_count
        full_accept=sum(x['method']=='full_smooth' for x in numeric)
        fallback_count=strict_count-full_accept
        need(len(evaluated_nodes)==counts['boxes_completed'],'closed evaluated-node accounting')
        need(counts['full_smooth_bundle_calls']==strict_count and counts['split_fallback_calls']==fallback_count and counts['split_bundle_calls']==boundary_count+fallback_count,'closed strict/boundary/fallback calls from all node gates')
    return numeric,{'full_tree_structure_checked_this_run':True,'tree_nodes':len(nodes),
                    'numeric_accepted_leaves':len(numeric),'unresolved_terminals':len(queue),
                    'seed_splits':seed,'midpoint_splits':mid,
                    'exact_cover_including_unresolved':True,'full_accepted_cover':closed,
                    'accepted_parameter_area':str(accepted_area),'unresolved_parameter_area':str(pending_area)}


def main():
    started=monotonic()
    parser=argparse.ArgumentParser()
    inputs=parser.add_mutually_exclusive_group()
    inputs.add_argument('--summary',type=Path)
    inputs.add_argument('--full',type=Path)
    parser.add_argument('--dependency-map',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    source=args.full or args.summary or ROOT/'certificate_summary.json'
    data=json.loads(source.read_text())
    mapping=json.loads(args.dependency_map.read_text()) if args.dependency_map else {}
    need(data['round']==18 and data['candidate_attempts']==1,'one round18 attempt')
    need(data['parameters']=={'R':'15/8','e_by_power':{str(p):str(e) for p,e in E.items()}},'unique candidate')
    need(data['status'] in ('not_closed','closed_pending_independent_review'),'reviewable run status')
    need(data['limits']=={'seconds_including_preparation_aux_cap_and_boxes':285.0,'boxes_started':30000,'width_floor':str(Q(1,2**60)),'tail_r_n_max':28,'tail_even_k_max':14},'fixed authorized limits and truncations')
    need(0<=data['timing']['preparation_aux_and_cap_seconds']<=data['timing']['inclusive_certification_seconds']<=285,'recorded inclusive certification time')
    minimal=data.get('certificate_format')=='minimal_hybrid_v1'
    need(sha256((ROOT/'hybrid_certify.py').read_bytes()).hexdigest()==data['certifier_sha256'],'certifier source hash')
    for name,wanted in data['dependencies'].items():resolve(name,wanted,mapping)
    c,base,polynomials,tail=prepare_primitives(data,mapping)
    auxiliary=auxiliary_and_cap(c,data,mapping)
    if args.full:
        numeric,structure=tree_check(data,minimal)
        if data['status']=='closed_pending_independent_review':
            need(structure['full_accepted_cover'],'closed attempt must cover root')
            selected=numeric;all_numeric=True
        else:
            indices=sorted({0,len(numeric)//2,len(numeric)-1}) if numeric else []
            selected=[numeric[i] for i in indices];all_numeric=False
    else:
        need(data['status']=='not_closed','compact summary cannot verify all leaves of a closed attempt')
        need('partition_tree' not in data and 'accepted_leaf_data' not in data,'summary mode does not claim a full tree')
        selected=data.get('selected_numeric_leaves',data.get('selected_accepted_leaves',[]))
        selected=[x for x in selected if x['method']!='analytic_cap']
        need(len(selected)<=3,'failed summary replay limited to three numeric leaves')
        structure={'full_tree_structure_checked_this_run':False,'reason':'Compact summary only; --full FILE is required for full tree coverage replay.'}
        all_numeric=False
    checks=[]
    for index,leaf in enumerate(selected,1):
        checks.append(numeric_leaf(c,base,polynomials,tail,leaf,minimal))
        if index%200==0:
            print(json.dumps({'verified_numeric_leaves':index,'total_selected':len(selected),'elapsed_seconds':monotonic()-started}),flush=True)
    result={'status':('closed_all_numeric_leaves_and_cap_independently_verified' if all_numeric else
                      'failed_attempt_full_tree_and_sampled_leaves_verified' if args.full else 'failed_attempt_summary_sampled_leaves_and_cap_verified'),
            'mode':'full' if args.full else 'summary','candidate_status':data['status'],
            'certificate_format':'minimal_hybrid_v1' if minimal else 'detailed_original',
            'verification_elapsed_seconds':monotonic()-started,
            **structure,**auxiliary,'numeric_leaves_independently_recomputed':len(checks),
            'all_numeric_leaf_formulas_recomputed':all_numeric,'checks':checks,
            'shared_reviewed_primitives':'round008 Fraction/log; round017 Poly, AD, shifted_bundles and tail_bundles',
            'certifier_bound_or_adaptive_functions_called':False,'new_adaptive_runs':0,'new_LPs':0,
            'input_sha256':sha256(source.read_bytes()).hexdigest(),
            'verifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Only a closed full-tree replay with every numeric leaf plus cap/auxiliary/analytic transfer audit may update a rigorous upper bound.'}
    output=args.output or ROOT/('hybrid_full_verification.json' if args.full else 'hybrid_summary_verification.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2),flush=True)


if __name__=='__main__':
    main()
