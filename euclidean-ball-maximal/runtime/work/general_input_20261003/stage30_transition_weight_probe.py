"""Actual matched-exit transition weights, full original masks.
Rational polynomial primitives with certified quadratic-root enclosures.
No spatial grid, no source subsampling, no fixed-H substitution.
"""
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
from bisect import bisect_right
from math import isqrt
import json,time,sys,hashlib
sys.dont_write_bytecode=True
from stage25_stopped_energy_probe import MomentMeasure
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'output/general_input_20261003'
BITS=96
ZERO=(F(0),F(0))
def plus(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def minus(a,b):return plus(a,neg(b))
def times(a,b):
    vals=[x*y for x in a for y in b];return min(vals),max(vals)
def scalar(a,c):return times(a,(c,c))
def enc(a):return dict(lo=str(a[0]),hi=str(a[1]),display=float((a[0]+a[1])/2),width=str(a[1]-a[0]))
def peval(p,x):
    out=F(0)
    for c in reversed(p):out=out*x+c
    return out
def ieval(p,x):
    out=ZERO
    for c in reversed(p):out=plus(times(out,x),(c,c))
    return out
def pmul(p,q):
    ans=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):ans[i+j]+=a*b
    return ans
def pint(p):return [F(0)]+[v/(i+1) for i,v in enumerate(p)]
def integ(p,l,h):return minus(ieval(pint(p),h),ieval(pint(p),l))
def sqrtbox(q,bits=BITS):
    assert q>=0
    sn=isqrt(q.numerator);sd=isqrt(q.denominator)
    if sn*sn==q.numerator and sd*sd==q.denominator:v=F(sn,sd);return v,v
    factor=2**bits;n=isqrt(q.numerator*factor*factor//q.denominator)
    lo=F(n,factor);hi=F(n+1,factor);assert lo*lo<q<hi*hi
    return lo,hi

def roots(p,t,L):
    a,b,c=p[2],p[1],p[0]-t
    if a==0:
        if b==0:return []
        z=-c/b;return [(z,z)] if 0<z<L else []
    disc=b*b-4*a*c
    if disc<0:return []
    if disc==0:
        z=-b/(2*a);return [(z,z)] if 0<z<L else []
    boxes=[]
    for sign in (-1,1):
        for bits in (BITS,192,384,768):
            s=sqrtbox(disc,bits);box=scalar(plus((-b,-b),scalar(s,F(sign))),1/(2*a))
            if box[1]<=0 or box[0]>=L:break
            if 0<box[0]<=box[1]<L:boxes.append(box);break
        else:raise AssertionError('Root overlaps piece boundary at 768bits')
    boxes.sort();assert all(boxes[i][1]<boxes[i+1][0] for i in range(len(boxes)-1))
    return boxes

def regions(p,levels,L):
    interior=[]
    for t in levels:interior.extend(roots(p,t,L))
    interior=sorted(set(interior))
    assert all(interior[i][1]<interior[i+1][0] for i in range(len(interior)-1))
    bounds=[(F(0),F(0))]+interior+[(L,L)]
    return [(l,h,(l[1]+h[0])/2) for l,h in zip(bounds,bounds[1:])],interior

class Linear:
    def __init__(self,events):
        self.events={x:v for x,v in events.items() if v};self.knots=sorted(self.events)
        self.s=[];self.b=[];s=b=F(0)
        for x in self.knots:s+=self.events[x];b+=x*self.events[x];self.s.append(s);self.b.append(b)
        assert s==b==0
    def value(self,x):
        i=bisect_right(self.knots,x)-1
        return F(0) if i<0 else self.s[i]*x-self.b[i]
    def slope(self,x):
        i=bisect_right(self.knots,x)-1
        return F(0) if i<0 else self.s[i]
    def save(self):return [dict(x=str(x),value=str(self.value(x)),slope_jump=str(self.events[x])) for x in self.knots]

def build_h(cells,a,kind):
    events=defaultdict(lambda:defaultdict(F))
    for row in cells:
        if kind=='eligible' and row['J']==1:continue
        if kind=='J1' and row['J']!=1:continue
        k=row['K'];r=F(row['radius']);l,h=F(row['lo']),F(row['hi'])
        factor=1/(2*r)
        for x,v in [(l-r,1),(l+r,-1),(h-r,-1),(h+r,1)]:events[k][x]+=v*factor
    return {k:Linear(e) for k,e in events.items()}

def audit(source,saved):
    atoms=[(F(v['location']),F(v['mass'])) for v in source]
    P=MomentMeasure(atoms=atoms)
    stop=saved['complete_actual_observer_stop'];alpha=F(saved['alpha']);a=F(stop['entrance_a']);D=stop['D']
    eta=MomentMeasure(cells=[tuple(F(v[k]) for k in ('lo','hi','density')) for v in stop['complete_eta_density']])
    u=lambda x:eta.left_distance(x)-P.left_distance(x)
    du=lambda x:eta.prefix(x)[0]-P.prefix(x)[0]
    measures={y for y,w in atoms}|{z for l,h,d in eta.cells for z in (l,h)}
    radii={k:F(1,2**k)/(2*a) for k in range(1,D+1)}
    ts={k:alpha*r*r/16 for k,r in radii.items()}
    outputs=[]
    for kind in ('full','eligible','J1'):
        hs=build_h(saved['complete_original_cells'],a,kind)
        target_by_k=defaultdict(F)
        for row in saved['complete_original_cells']:
            if kind=='eligible' and row['J']==1:continue
            if kind=='J1' and row['J']!=1:continue
            target_by_k[row['K']]+=F(row['matched_signed_residual_spatial_integral'])
        layers=[];total_source=total_original=F(0);total_eta=total_energy=ZERO;root_count=piece_count=0
        for k,hfun in sorted(hs.items()):
            t=ts[k]
            original=sum((w*hfun.value(y) for y,w in atoms),F(0))
            lost=sum((w*hfun.value(y)*min(u(y)/t,F(1)) for y,w in atoms),F(0))
            assert 0<=lost<=original
            cuts=sorted(measures|set(hfun.knots));out_eta=out_energy=ZERO;parts=[]
            for lo,hi in zip(cuts,cuts[1:]):
                mid=(lo+hi)/2;hc=hfun.value(lo);slope=hfun.slope(mid)
                if hc==slope==0:continue
                idx=bisect_right(eta.points,mid)-1
                density=eta.cells[idx][2] if idx>=0 and mid<eta.cells[idx][1] else F(0)
                poly=[u(lo),du(lo),density/2];L=hi-lo
                assert peval(poly,L)==u(hi) and min(poly[0],peval(poly,L))>=0
                minimum=min(poly[0],peval(poly,L))
                if poly[2] and 0<-poly[1]/(2*poly[2])<L:
                    minimum=min(minimum,poly[0]-poly[1]*poly[1]/(4*poly[2]))
                hpoly=[hc,slope];dpoly=[poly[1],2*poly[2]]
                sub,rr=regions(poly,[t],L);root_count+=len(rr)
                ei=zi=ZERO
                for left,right,center in sub:
                    below=peval(poly,center)<t
                    epoly=[v*density/t for v in pmul(hpoly,poly)] if below else [v*density for v in hpoly]
                    zpoly=[v/t for v in pmul(hpoly,pmul(dpoly,dpoly))] if below else [F(0)]
                    ei=plus(ei,integ(epoly,left,right));zi=plus(zi,integ(zpoly,left,right))
                # Both true integrands are nonnegative. Clamp harmless interval roundoff.
                ei=(max(F(0),ei[0]),ei[1]);zi=(max(F(0),zi[0]),zi[1])
                assert 0<=ei[0]<=ei[1] and 0<=zi[0]<=zi[1]
                out_eta=plus(out_eta,ei);out_energy=plus(out_energy,zi);piece_count+=1
                parts.append(dict(lo=str(lo),hi=str(hi),u_coefficients=list(map(str,poly)),h_coefficients=list(map(str,hpoly)),
                    eta_density=str(density),exact_minimum_u_on_piece=str(minimum),
                    whole_piece_saturated=minimum>=t,threshold_roots_from_lo=[enc(r) for r in rr],
                    eta_transition_integral=enc(ei),energy_integral=enc(zi)))
            got=minus((lost,lost),plus(out_eta,out_energy));target=target_by_k[k]
            assert got[0]<=target<=got[1] and got[1]-got[0]<F(1,10**12)
            total_original+=original;total_source+=lost;total_eta=plus(total_eta,out_eta);total_energy=plus(total_energy,out_energy)
            layers.append(dict(k=k,radius=str(radii[k]),height=str(t),h_knots=hfun.save(),original_source_pairing=str(original),
                source_retained=str(original-lost),source_transition_loss=str(lost),eta_transition_pairing=enc(out_eta),energy_pairing=enc(out_energy),
                source_minus_eta_minus_energy=enc(got),stage29_matched_residual_exact=str(target),agreement_certified=True,space_pieces=parts))
        got=minus((total_source,total_source),plus(total_eta,total_energy));target=sum(target_by_k.values(),F(0))
        assert got[0]<=target<=got[1] and got[1]-got[0]<F(1,10**11)
        # Complete transition overlap: linear Q_j restricted by quadratic u.
        overlaps=[];weighted=F(0);global_max=ZERO
        for j in range(1,D):
            qe=defaultdict(F)
            for k,hfun in hs.items():
                if k>j:
                    for x,v in hfun.events.items():qe[x]+=v
            Q=Linear(qe)
            if not Q.knots:continue
            cuts=sorted(measures|set(Q.knots));candidates=[];accepted=0
            for lo,hi in zip(cuts,cuts[1:]):
                mid=(lo+hi)/2;hc=Q.value(lo);slope=Q.slope(mid)
                if hc==slope==0:continue
                idx=bisect_right(eta.points,mid)-1
                density=eta.cells[idx][2] if idx>=0 and mid<eta.cells[idx][1] else F(0)
                poly=[u(lo),du(lo),density/2];L=hi-lo
                sub,rr=regions(poly,[ts[j+1],ts[j]],L)
                for left,right,center in sub:
                    if not ts[j+1]<peval(poly,center)<ts[j]:continue
                    accepted+=1
                    for box in (left,right):
                        value=ieval([hc,slope],box)
                        candidates.append(dict(position=enc(plus((lo,lo),box)),value=enc(value)))
            if candidates:
                ml=max(F(v['value']['lo']) for v in candidates);mh=max(F(v['value']['hi']) for v in candidates)
                best=(max(F(0),ml),max(F(0),mh))
            else:best=ZERO
            contribution=sum((w*Q.value(y) for y,w in atoms if ts[j+1]<u(y)<ts[j]),F(0))
            weighted+=contribution;global_max=(max(global_max[0],best[0]),max(global_max[1],best[1]))
            overlaps.append(dict(j=j,height_lower=str(ts[j+1]),height_upper=str(ts[j]),
                supremum_Q_on_open_transition=enc(best),nonempty_transition_pieces=accepted,
                source_weighted_transition_exact=str(contribution),all_supremum_candidates=candidates))
        outputs.append(dict(observer_part=kind,original_source_pairing=str(total_original),source_retained=str(total_original-total_source),
            source_transition_loss=str(total_source),eta_transition_pairing=enc(total_eta),energy_pairing=enc(total_energy),
            source_minus_eta_minus_energy=enc(got),stage29_residual_exact=str(target),agreement_certified=True,
            threshold_root_count=root_count,space_piece_count=piece_count,layers=layers,
            all_h_support_pieces_saturated=all(c['whole_piece_saturated'] for layer in layers for c in layer['space_pieces']),
            transition_overlap=dict(global_supremum_Q=enc(global_max),source_weighted_integral_exact=str(weighted),per_transition=overlaps)))
    assert sum(F(o['stage29_residual_exact']) for o in outputs[1:])==F(outputs[0]['stage29_residual_exact'])
    assert F(outputs[0]['stage29_residual_exact'])==F(saved['full_spatial_integrals']['matched_signed_residual'])
    return dict(label=saved['label'],scope=saved['scope'],b_over_a=saved['b_over_a'],alpha=saved['alpha'],
        original_volume=saved['original_full_volume'],eligible_volume=saved['eligible_volume'],J1_volume=saved['J1_volume'],
        actual_D_cutoff=D,complete_original_observer_cells=saved['complete_original_cells'],actual_eta=stop['complete_eta_density'],
        original_H_segments=stop['H_linear_segments'],parts=outputs)

def main():
    start=time.monotonic();raw=(DATA/'stage29_matched_height_probe.json').read_bytes();base=json.loads(raw)
    rounds=[]
    for tag in ('frozen_ratio4_rounds','new_ratio8_rounds'):
        for rd in base[tag]:
            scopes=[]
            for saved in rd['scopes']:
                got=audit(rd['source'],saved);scopes.append(got);o=got['parts'][0]
                print(tag,rd['round'],got['scope'],'S',o['stage29_residual_exact'],'energy',o['energy_pairing']['display'],
                    'maxQ',o['transition_overlap']['global_supremum_Q']['display'],
                    'sourceQ',float(F(o['transition_overlap']['source_weighted_integral_exact'])),flush=True)
            rounds.append(dict(input_group=tag,round=rd['round'],source=rd['source'],scopes=scopes))
    out=dict(status='passed',dimension=1,source_file='stage29_matched_height_probe.json',source_sha256=hashlib.sha256(raw).hexdigest(),
        method='Fraction polynomial primitives; quadratic threshold roots enclosed by exact integer-square-root dyadic intervals; interval Horner evaluates every primitive and transition maximum.',
        nominal_root_bits=BITS,rounds=rounds,no_observer_source_or_level_grid_sampling=True,
        scope_count=sum(len(r['scopes']) for r in rounds),
        h_support_piece_count_including_separate_ledgers=sum(part['space_piece_count'] for r in rounds for scope in r['scopes'] for part in scope['parts']),
        energy_threshold_root_count=sum(part['threshold_root_count'] for r in rounds for scope in r['scopes'] for part in scope['parts']),
        all_supremum_candidate_count=sum(len(z['all_supremum_candidates']) for r in rounds for scope in r['scopes'] for part in scope['parts'] for z in part['transition_overlap']['per_transition']),
        full_original_J1_and_eligible_saved_separately=True,main_weak_theorem_proved=False,
        finite_tests_do_not_prove_general_transition_budget=True,elapsed_seconds=time.monotonic()-start)
    (DATA/'stage30_transition_weight_probe.json').write_text(json.dumps(out,indent=2)+'\n')
    print('completed',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
