#!/usr/bin/env python3
"""Exact finite-n1 common-threshold flow and record-kernel stress tests.
Writes only requested outputs, imports archived046 constructors dynamically.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from bisect import bisect_left,bisect_right
import importlib.util,argparse,json,random,time,hashlib,sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[2]
KERNEL_LABELS={'minimal_four_atoms','four_atom_perturbation_s471016','geometric_contact_clouds_s491101','geometric_contact_clouds_s491107','geometric_contact_clouds_s491108','chain_L2_alpha1','chain_L4_alpha1','chain_L8_alpha1','chain_L4_alpha2','chain_L8_alpha2','chain_L16_alpha2','chain_L4_alpha3/2','chain_L8_alpha3/2'}

def load_module(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def cases(old,prefix):
    minimal=[(F(0),F(13,256)),(F(1,16),F(3,64)),(F(3,32),F(9,128)),(F(10),F(213,256))]
    yield 0,'minimal_four_atoms',minimal,F(1),dict(family='minimal_counterexample')
    for i in range(11):
        rng=random.Random(471007+i);weights=[F(rng.randint(5,26),256) for _ in range(3)]
        positions=[F(0),F(1,16)+F(rng.randint(-2,2),512),F(3,32)+F(rng.randint(-2,2),512)]
        atoms=sorted(list(zip(positions,weights))+[(F(10),1-sum(weights))]);alpha=rng.choice((F(3,4),F(1),F(3,2),F(2)))
        yield 0,f'four_atom_perturbation_s{471007+i}',atoms,alpha,dict(seed=471007+i,family='dyadic_four_atom_perturbation')
    counts={1:0,2:0}
    for b,label,atoms,alpha,params in prefix.cases(471100):
        if b in counts and counts[b]<8 and len(atoms)<=12:
            counts[b]+=1;yield 1,label,atoms,alpha,dict(params,search_family=b)
    assert counts=={1:8,2:8}
    for alpha in (F(1),F(3,2),F(2)):
        for L in (2,4,8,16,32):
            yield 2,f'chain_L{L}_alpha{alpha}',old.old.old.chain(L,den=4,position=F(9,8),weight=F(27,16)),alpha,dict(family='chain_depth',L=L,alpha=str(alpha))

def context(old,atoms,alpha):
    saved,_,radii,_=old.observer(atoms,alpha,'band',4);D=saved['D'];lo=alpha/4;hi=alpha/2
    edges=sorted({y+sign*radii[k] for y,w in atoms for k in range(1,D+1) for sign in (-1,1)})
    cells=[];cuts={lo,hi}
    for c in saved['observer_cells']:
        a=F(c['lo']);b=F(c['hi']);nodes=[a]+edges[bisect_right(edges,a):bisect_left(edges,b)]+[b]
        for l,r in zip(nodes,nodes[1:]):
            mid=(l+r)/2;trace=[sum((w for y,w in atoms if abs(y-mid)<radii[k]),F(0))/(2*radii[k]) for k in range(1,D+1)]
            assert trace[0]<=alpha/8 and max(trace)>hi
            cells.append((l,r,c['J'],trace));cuts.update(g for g in trace if lo<g<hi)
    return saved,radii,cells,sorted(cuts)

def linear_sign_integrals(v0,v1,width):
    net=width*(v0+v1)/2
    if min(v0,v1)>=0:return net,F(0)
    if max(v0,v1)<=0:return F(0),-net
    pos=width*max(v0,v1)**2/(2*abs(v1-v0));return pos,pos-net

def flow_probe(old,record,batch,label,atoms,alpha,params):
    began=time.monotonic();saved,radii,cells,cuts=context(old,atoms,alpha)
    X=alpha*F(saved['eligible_volume']);assert X>0
    totals={k:F(0) for k in ('M','Q','diag','B','F','flow_pos','flow_neg')};beta_positive=beta_negative=F(0)
    lo=alpha/4;hi=alpha/2;first_positive=None;first_negative=None;pair_entries=0;edge_piece_count=0;min_beta_flow=None;max_beta_flow=None
    for a,b in zip(cuts,cuts[1:]):
        beta=(a+b)/2;assigned=[]
        for l,r,J,trace in cells:
            k=next(k for k,g in enumerate(trace,1) if g>beta);g=trace[k-1]
            assert J<=k and beta<g<=2*beta<=alpha
            assigned.append(dict(lo=str(l),hi=str(r),J=J,K=k,radius=str(radii[k])))
        hs=old.build_h(assigned,alpha/8,'eligible');A=old.add(hs.values())
        sourceA={y:A.value(y) for y,w in atoms}
        M=sum((w*sourceA[y] for y,w in atoms),F(0));Q=sum((w*sourceA[y]**2 for y,w in atoms),F(0));diag=sum((w*sum((h.value(y)**2 for h in hs.values()),F(0)) for y,w in atoms),F(0))
        prefixes={k:old.add(hs[j] for j in hs if j<k) for k in hs};B=pos=neg=F(0)
        for row,cell in zip(assigned,cells):
            l,r,J,trace=cell;k=row['K'];sk=prefixes[k];rad=radii[k]
            nodes=[l]+sk.knots[bisect_right(sk.knots,l):bisect_left(sk.knots,r)]+[r]
            B+=trace[k-1]*sum(((v-u)*(sk.value(u)+sk.value(v))/2 for u,v in zip(nodes,nodes[1:])),F(0))
            for y,w in atoms:
                u=max(l,y-rad);v=min(r,y+rad)
                if u>=v:continue
                sy=sk.value(y);nodes=[u]+sk.knots[bisect_right(sk.knots,u):bisect_left(sk.knots,v)]+[v]
                for s,t in zip(nodes,nodes[1:]):
                    p,n=linear_sign_integrals(sy-sk.value(s),sy-sk.value(t),t-s);pos+=p*w/(2*rad);neg+=n*w/(2*rad);edge_piece_count+=1
        flow=(Q-diag)/2-B;assert flow==pos-neg and diag<=M and beta*X/alpha<=M<=2*beta*X/alpha
        if min_beta_flow is None or flow<min_beta_flow:min_beta_flow=flow
        if max_beta_flow is None or flow>max_beta_flow:max_beta_flow=flow
        assert B<=2*beta*X/alpha
        weight=(b-a)/(hi-lo)
        for key,value in dict(M=M,Q=Q,diag=diag,B=B,F=flow,flow_pos=pos,flow_neg=neg).items():totals[key]+=weight*value
        if flow>0:
            beta_positive+=b-a
            if first_positive is None:first_positive=dict(lo=str(a),hi=str(b),F=str(flow),F_over_X=str(flow/X))
        if flow<0:
            beta_negative+=b-a
            if first_negative is None:first_negative=dict(lo=str(a),hi=str(b),F=str(flow),F_over_X=str(flow/X))
    selected=[cells[i][3] for i in sorted({i*(len(cells)-1)//5 for i in range(6)})]
    for x in selected:
        for z in selected:pair_entries+=record.pair_audit(x,z,lo,hi)
    assert totals['F']==totals['flow_pos']-totals['flow_neg']
    assert F(7,48)*X*X<=totals['Q']<=F(9,4)*X+2*totals['F']
    return dict(batch=batch,label=label,params=params,N=len(atoms),alpha=str(alpha),D=saved['D'],X=str(X),observer_cell_count=len(cells),beta_interval_count=len(cuts)-1,mean={key:str(v) for key,v in totals.items()},ratio={key+'_over_X':str(v/X) for key,v in totals.items()},beta_positive_flow_probability=str(beta_positive/(hi-lo)),beta_negative_flow_probability=str(beta_negative/(hi-lo)),flow_cancellation_fraction=str(2*min(totals['flow_pos'],totals['flow_neg'])/(totals['flow_pos']+totals['flow_neg'])) if totals['flow_pos']+totals['flow_neg'] else '0',first_positive_beta_interval=first_positive,first_negative_beta_interval=first_negative,minimum_beta_F_over_X=str(min_beta_flow/X),maximum_beta_F_over_X=str(max_beta_flow/X),pair_entries_checked=pair_entries,edge_affine_piece_count=edge_piece_count,atoms=[dict(location=str(y),mass=str(w)) for y,w in atoms],exact_error='0',elapsed_seconds=time.monotonic()-began)

def kernel_probe(old,record,atoms,alpha):
    """Full E² positive/negative L integration; x,z cells fix all ball members.

    Inside |x-z|<Rj the kernel is nonpositive; outside it is nonnegative.
    W/|I| weights one common beta. 'record_overlap_removed' retains nonempty
    individual record labels but removes their common-beta overlap entirely.
    """
    began=time.monotonic();saved,radii,cells,cuts=context(old,atoms,alpha);lo=alpha/4;hi=alpha/2
    data=[]
    for xa,xb,J,trace in cells:
        intervals=record.intervals(trace,lo,hi);mid=(xa+xb)/2
        members={k:frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-mid)<radii[k]) for k in range(1,saved['D']+1)}
        data.append((xa,xb,trace,intervals,members))
    pos=neg=unweighted=charged=exclusive_charged=F(0);checked=weighted_pairs=0;pairwise_positive_net=pairwise_negative_net=F(0);witness=None;annular_pairs=set()
    for xa,xb,tx,ix,mx in data:
        for za,zb,tz,iz,mz in data:
            rectangle=(xb-xa)*(zb-za)
            for j,(la,lb) in ix.items():
                R=radii[j]
                square=lambda t:max(F(0),t)**2
                strip=sum((sign*(square(zb+shift)-square(za+shift))/2 for shift,sign in ((R-xa,1),(R-xb,-1),(-R-xa,-1),(-R-xb,1))),F(0))
                assert 0<=strip<=rectangle
                for k,(ra,rb) in iz.items():
                    if j>=k:continue
                    checked+=1;r=radii[k]
                    first=sum((atoms[i][1] for i in mx[j]&mz[k]),F(0))/(4*R*r)
                    second=tz[k-1]/(2*R);assert 0<=first<=second
                    p=first*(rectangle-strip);n=(second-first)*strip
                    overlap=max(F(0),min(lb,rb)-max(la,ra))
                    densitygap=max(F(0),tx[j-1]-tz[j-1]);exclusive=mx[j]-mz[j]
                    exclusive_density=sum((atoms[i][1] for i in exclusive),F(0))/(2*R)
                    assert overlap<=densitygap<=exclusive_density
                    unweighted+=p;charged+=densitygap*p/(hi-lo)
                    exclusive_charged+=exclusive_density*p/(hi-lo)
                    if p:
                        for yi in mx[j]&mz[k]:
                            for vi in exclusive:
                                distance=abs(atoms[yi][0]-atoms[vi][0]);assert R-r<distance<2*R
                                annular_pairs.add((j,k,yi,vi))
                    if overlap:
                        weighted_pairs+=1;factor=overlap/(hi-lo);pos+=factor*p;neg+=factor*n
                        if p>n:pairwise_positive_net+=factor*(p-n)
                        if n>p:pairwise_negative_net+=factor*(n-p)
                        if p>n and witness is None:witness=dict(xlo=str(xa),xhi=str(xb),zlo=str(za),zhi=str(zb),j=j,k=k,W=str(overlap),density_difference=str(densitygap),first_density=str(first),second_density=str(second),rectangle_area=str(rectangle),coarse_strip_area=str(strip),weighted_positive_cost=str(factor*p),weighted_negative_cost=str(factor*n))
    assert pos<=unweighted and pos<=charged<=exclusive_charged and pos-neg==pairwise_positive_net-pairwise_negative_net
    X=alpha*F(saved['eligible_volume'])
    return dict(status='passed',domain='entire E_D x E_D',common_beta_weight='W_jk / |I|',record_overlap_removed_scope='j<k with both individual record intervals nonempty; replace W/|I| by1',positive=str(pos),negative=str(neg),net=str(pos-neg),positive_over_X=str(pos/X),negative_over_X=str(neg/X),record_overlap_removed_positive=str(unweighted),record_overlap_removed_positive_over_X=str(unweighted/X),density_difference_charged_positive=str(charged),density_difference_charged_positive_over_X=str(charged/X),exclusive_second_source_charged_positive=str(exclusive_charged),exclusive_second_source_charged_positive_over_X=str(exclusive_charged/X),exclusive_second_source_annular_support='R_j-R_k < |y-v| < 2 R_j',exclusive_second_source_annular_pairs_checked=len(annular_pairs),kernel_cancellation_fraction=str(2*min(pos,neg)/(pos+neg)) if pos+neg else '0',checked_rectangle_label_pairs=checked,weighted_rectangle_label_pairs=weighted_pairs,pairwise_positive_net=str(pairwise_positive_net),pairwise_negative_net=str(pairwise_negative_net),positive_rectangle_witness=witness,exact_error='0',elapsed_seconds=time.monotonic()-began)

def record_energy_probe(old,record,atoms,alpha):
    """Independent full E² record-collapse energies and symmetric kernel.

    Ef/Er use trace running maxima directly. Their width/square expressions
    are independently compared with individual first-exit record intersections.
    Gsym integrates H-(uj*gk+uk*gj)/2 on source rectangles, never from F.
    """
    began=time.monotonic();saved,radii,cells,_=context(old,atoms,alpha);lo=alpha/4;hi=alpha/2;I=hi-lo
    data=[]
    for xa,xb,J,trace in cells:
        ix=record.intervals(trace,lo,hi);mid=(xa+xb)/2;running=[F(0)]
        for g in trace:running.append(max(running[-1],g))
        clipped=[min(hi,max(lo,g)) for g in running]
        assert sum((clipped[k]**2-clipped[k-1]**2 for k in range(1,len(clipped))),F(0))==hi**2-lo**2
        members={k:frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-mid)<radii[k]) for k in ix}
        data.append((xa,xb,trace,ix,running,members))
    Ef=Er=B=Tobs=H=Gsym=F(0);directB=directT=F(0);forward_checks=reverse_checks=pair_checks=0
    for xa,xb,tx,ix,rx,mx in data:
        for za,zb,tz,iz,rz,mz in data:
            rectangle=(xb-xa)*(zb-za);cache={}
            def strip(k):
                if k not in cache:
                    R=radii[k];square=lambda t:max(F(0),t)**2
                    area=sum((sign*(square(zb+shift)-square(za+shift))/2 for shift,sign in ((R-xa,1),(R-xb,-1),(-R-xa,-1),(-R-xb,1))),F(0))
                    assert 0<=area<=rectangle;cache[k]=area
                return cache[k]
            # Collapse all later z labels at fixed x record j.
            for j,(a,b) in ix.items():
                lower=max(lo,rx[j-1],rz[j]);upper=min(hi,tx[j-1])
                square_weight=max(F(0),upper**2-lower**2);direct_square=width_cost=F(0)
                assert square_weight<=min(hi,max(lo,rx[j]))**2-min(hi,max(lo,rx[j-1]))**2
                for k,(c,d) in iz.items():
                    if j>=k:continue
                    u=max(a,c);v=min(b,d)
                    if u<v:direct_square+=v*v-u*u;width_cost+=(v-u)*tz[k-1]
                assert direct_square==square_weight and width_cost<=square_weight<=2*width_cost
                forward_checks+=1;factor=strip(j)/(2*radii[j]*I);Ef+=square_weight*factor;B+=width_cost*factor
            # Collapse all earlier x labels at fixed z record k.
            for k,(c,d) in iz.items():
                lower=max(lo,rz[k-1]);upper=min(hi,tz[k-1],rx[k-1])
                square_weight=max(F(0),upper**2-lower**2);direct_square=width_cost=F(0)
                assert square_weight<=min(hi,max(lo,rz[k]))**2-min(hi,max(lo,rz[k-1]))**2
                for j,(a,b) in ix.items():
                    if j>=k:continue
                    u=max(a,c);v=min(b,d)
                    if u<v:direct_square+=v*v-u*u;width_cost+=(v-u)*tx[j-1]
                assert direct_square==square_weight and width_cost<=square_weight<=2*width_cost
                reverse_checks+=1;factor=strip(k)/(2*radii[k]*I);Er+=square_weight*factor;Tobs+=width_cost*factor
            # Direct symmetric kernel integral, without using the flow identity.
            for j,(a,b) in ix.items():
                for k,(c,d) in iz.items():
                    if j>=k:continue
                    width=max(F(0),min(b,d)-max(a,c))
                    if not width:continue
                    pair_checks+=1;R=radii[j];r=radii[k]
                    first=sum((atoms[i][1] for i in mx[j]&mz[k]),F(0))/(4*R*r)
                    h_integral=first*rectangle
                    a_integral=tz[k-1]*strip(j)/(2*R)
                    b_integral=tx[j-1]*strip(k)/(2*r)
                    factor=width/I;H+=factor*h_integral;directB+=factor*a_integral;directT+=factor*b_integral
                    Gsym+=factor*(h_integral-(a_integral+b_integral)/2)
    X=alpha*F(saved['eligible_volume'])
    assert B==directB and Tobs==directT and B<=Ef<=2*B and Tobs<=Er<=2*Tobs
    assert Ef<=F(3,4)*X and Er<=F(3,4)*X and Ef+Er<=F(3,2)*X
    imbalance=(Tobs-B)/2;assert abs(imbalance)<=F(3,8)*X
    return dict(status='passed',domain='entire E_D x E_D',threshold_interval=[str(lo),str(hi)],Ef=str(Ef),Er=str(Er),Ef_over_X=str(Ef/X),Er_over_X=str(Er/X),energy_sum_over_X=str((Ef+Er)/X),Bbar=str(B),Tobsbar=str(Tobs),Bbar_over_X=str(B/X),Tobsbar_over_X=str(Tobs/X),joint_Hbar=str(H),Gsym=str(Gsym),Gsym_over_X=str(Gsym/X),imbalance=str(imbalance),imbalance_over_X=str(imbalance/X),symmetric_kernel_definition='H-(u_j(x-z)g_k(z)+u_k(x-z)g_j(x))/2; direct source rectangle integration',forward_record_collapse_checks=forward_checks,reverse_record_collapse_checks=reverse_checks,weighted_symmetric_kernel_pairs_checked=pair_checks,forward_bounds=True,reverse_bounds=True,energy_sum_bound=True,imbalance_bound=True,exact_error='0',elapsed_seconds=time.monotonic()-began)

def finalize(out):
    out['batches']=[]
    for b in sorted({c['batch'] for c in out['cases']}):
        group=[c for c in out['cases'] if c['batch']==b];maximum={}
        for key in ('F_over_X','flow_pos_over_X','flow_neg_over_X','Q_over_X'):
            c=max(group,key=lambda c:F(c['ratio'][key]));maximum[key]=dict(label=c['label'],exact=c['ratio'][key],display=float(F(c['ratio'][key])))
        minc=min(group,key=lambda c:F(c['ratio']['F_over_X']))
        out['batches'].append(dict(batch=b,count=len(group),positive_count=sum(F(c['mean']['F'])>0 for c in group),negative_count=sum(F(c['mean']['F'])<0 for c in group),zero_count=sum(F(c['mean']['F'])==0 for c in group),maximum=maximum,minimum_F_over_X=dict(label=minc['label'],exact=minc['ratio']['F_over_X'],display=float(F(minc['ratio']['F_over_X'])))))
    out['exact_checks']=dict(error='0',case_count=len(out['cases']),edge_affine_pieces=sum(c['edge_affine_piece_count'] for c in out['cases']),record_pair_entries=sum(c['pair_entries_checked'] for c in out['cases']),mean_F_equals_mean_positive_minus_negative=True,all_moment_bounds=True,full_Esquare_kernel_checks=sum('record_kernel' in c for c in out['cases']),record_energy_cases=sum('record_energy' in c for c in out['cases']))
    out['sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True);p.add_argument('--labels',nargs='*');args=p.parse_args()
    record=load_module(args.repo_root/'rounds/046/threshold_probe.py','record47');prefix=load_module(args.repo_root/'rounds/046/prefix_probe.py','prefix47');old=record.load(args.repo_root,args.output)
    start=time.monotonic();out=dict(status='running',dimension=1,arithmetic='Fraction exact',scope='actual original band with one common uniform threshold',cases=[],finite_tests_prove_general_averaged_bound=False)
    for batch,label,atoms,alpha,params in cases(old,prefix):
        if args.labels is not None and label not in args.labels:continue
        c=flow_probe(old,record,batch,label,atoms,alpha,params)
        if label in KERNEL_LABELS:
            c['record_kernel']=kernel_probe(old,record,atoms,alpha);assert c['record_kernel']['net']==c['mean']['F']
            energy=record_energy_probe(old,record,atoms,alpha)
            assert F(energy['Bbar'])==F(c['mean']['B'])
            assert F(c['mean']['F'])==F(energy['Gsym'])+F(energy['imbalance'])
            assert F(energy['joint_Hbar'])==(F(c['mean']['Q'])-F(c['mean']['diag']))/2
            c['record_energy']=energy
        out['cases'].append(c);args.output.write_text(json.dumps(out)+'\n')
        print(batch,label,'F/X',float(F(c['ratio']['F_over_X'])),'flow+/X',float(F(c['ratio']['flow_pos_over_X'])),'beta',c['beta_interval_count'],'sec',round(c['elapsed_seconds'],2),'kernel_sec',round(c.get('record_kernel',{}).get('elapsed_seconds',0),2),flush=True)
    finalize(out);out.update(status='passed',elapsed_seconds=time.monotonic()-start);args.output.write_text(json.dumps(out,indent=2)+'\n');print('PASSED',len(out['cases']),'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
