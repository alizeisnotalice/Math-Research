#!/usr/bin/env python3
"""Exact stopped dual-source testing, original scopes and recomputed low J/K.
Three frozen digital inputs. Common-c Brownian n1 entry stops, complete source
and radial integrals; H is an exact compact piecewise-linear signed function.
"""
from fractions import Fraction as F
from collections import defaultdict
from bisect import bisect_right
from pathlib import Path
import json,time,sys
sys.dont_write_bytecode=True
from stage18_stopped_flux_probe import full_E,difference,merge,flat,inside


class LinearH:
    def __init__(self,ex,en,R):
        changes=defaultdict(F)
        for sign,marks in [(1,ex),(-1,en)]:
            for j in range(1,len(R)):
                factor=F(sign)/(2*R[j])
                for l,h in zip(marks[j][::2],marks[j][1::2]):
                    for x,c in [(l-R[j],1),(l+R[j],-1),(h-R[j],-1),(h+R[j],1)]:
                        changes[x]+=c*factor
        self.knots=sorted(x for x,c in changes.items() if c)
        self.slopes=[];self.moments=[];self.values=[];s=b=F(0)
        for x in self.knots:
            s+=changes[x];b+=x*changes[x]
            self.slopes.append(s);self.moments.append(b);self.values.append(s*x-b)
        assert s==b==0
        self.segments=list(zip(self.knots,self.knots[1:],self.values,self.values[1:]))
        self.areas=[F(0)]
        for lo,hi,v,w in self.segments:self.areas.append(self.areas[-1]+(hi-lo)*(v+w)/2)
        assert sum(((hi-lo)*(v+w)/2 for lo,hi,v,w in self.segments),F(0))==0
    def __call__(self,x):
        i=bisect_right(self.knots,x)-1
        return F(0) if i<0 else self.slopes[i]*x-self.moments[i]
    def positive_integral(self,t=F(0)):
        out=F(0)
        for lo,hi,v,w in self.segments:
            v-=t;w-=t;length=hi-lo
            if min(v,w)>=0:out+=length*(v+w)/2
            elif max(v,w)<=0:continue
            elif v>0:out+=length*v*v/(2*(v-w))
            else:out+=length*w*w/(2*(w-v))
        return out
    def primitive(self,x):
        i=bisect_right(self.knots,x)-1
        if i<0:return F(0)
        if i>=len(self.knots)-1:return self.areas[-1]
        d=x-self.knots[i]
        return self.areas[i]+d*self.values[i]+self.slopes[i]*d*d/2
    def cap_upper(self,c):
        heights=sorted(set(v for v in self.values if v>0))
        candidates={F(0)}
        if heights:
            for i in range(16):candidates.add(heights[i*(len(heights)-1)//15])
        trials=[(t,t+c*self.positive_integral(t)) for t in sorted(candidates)]
        t,best=min(trials,key=lambda row:row[1])
        return best,t,trials


def observer(atoms,alpha,scope,ratio):
    b=alpha/2;a=b/ratio;geom=2*a
    E=full_E(atoms,alpha)
    if scope=='band':E=difference(E,full_E(atoms,2*alpha))
    D=1
    while F(1,2**D)>min(w for y,w in atoms)/(2*ratio):D+=1
    p=[F(1,2**j) for j in range(D+1)];R=[pp/geom for pp in p]
    events=defaultdict(list)
    for j,r in enumerate(R):
        for y,w in atoms:events[y-r].append((j,w));events[y+r].append((j,-w))
    cuts=sorted(set(events)|set(flat(E)));m=[F(0)]*(D+1);endpoints=flat(E)
    ex=defaultdict(list);en=defaultdict(list);Fmarks=defaultdict(list);cells=[];volume=j1=V=T=F(0)
    for lo,hi in zip(cuts,cuts[1:]):
        for j,w in events[lo]:m[j]+=w
        if not inside((lo+hi)/2,endpoints):continue
        J=next(j for j in range(1,D+1) if m[j]>p[j])
        K=next(j for j in range(1,D+1) if m[j]>ratio*p[j]);width=hi-lo;volume+=width
        assert J<=K
        if J==1:j1+=width;continue
        V+=width;gK=a*m[K]/p[K];gprev=a*m[J-1]/p[J-1]
        assert b<gK<=2*b and a/2<gprev<=a
        T+=width*(gK-gprev);ex[K].append((lo,hi));en[J-1].append((lo,hi))
        for j in range(J,K+1):Fmarks[j].append((lo,hi))
        cells.append(dict(lo=str(lo),hi=str(hi),J=J,K=K))
    assert volume==sum((hi-lo for lo,hi in E),F(0))==j1+V
    assert (b-a)*V<=T<=2*b*V
    ex=[flat(merge(ex[j])) for j in range(D+1)];en=[flat(merge(en[j])) for j in range(D+1)]
    H=LinearH(ex,en,R)
    sourceH=sum((w*H(y) for y,w in atoms),F(0));assert sourceH==T
    Hplus=H.positive_integral();assert Hplus<=V
    # Independently construct H from every actual F layer and compare all knots.
    terms=defaultdict(F)
    for j in range(2,D+1):
        for lo,hi in merge(Fmarks[j]):
            for sign,r in [(1,R[j]),(-1,R[j-1])]:
                for x,c in [(lo-r,1),(lo+r,-1),(hi-r,-1),(hi+r,1)]:terms[x]+=sign*c/(2*r)
    other=defaultdict(F)
    for sign,marks in [(1,ex),(-1,en)]:
        for j in range(1,D+1):
            for lo,hi in zip(marks[j][::2],marks[j][1::2]):
                for x,c in [(lo-R[j],1),(lo+R[j],-1),(hi-R[j],-1),(hi+R[j],1)]:other[x]+=sign*c/(2*R[j])
    assert {x:c for x,c in terms.items() if c}=={x:c for x,c in other.items() if c}
    return dict(scope=scope,alpha=str(alpha),entrance_a=str(a),terminal_b=str(b),b_over_a=ratio,D=D,
      original_intervals=[[str(lo),str(hi)] for lo,hi in E],original_volume=str(volume),alphaV=str(alpha*volume),
      J1_volume=str(j1),eligible_volume=str(V),T_exact=str(T),observer_cells=cells,
      all_original_observers_have_K_at_most_D=True,H_source_integral_exact=str(sourceH),
      H_positive_integral_exact=str(Hplus),H_negative_integral_exact=str(Hplus),
      H_min_exact=str(min(H.values,default=F(0))),H_max_exact=str(max(H.values,default=F(0))),
      H_knots=[dict(x=str(x),H=str(v)) for x,v in zip(H.knots,H.values)],
      H_linear_segments=[dict(lo=str(lo),hi=str(hi),Hlo=str(v),Hhi=str(w),slope=str((w-v)/(hi-lo)))
        for lo,hi,v,w in H.segments],
      H_both_actual_mask_formulas_match=True),en,R,H


def stopped_source(y,en,R,H):
    # Only actual entrance-mark layers and final1 can stop. Unmarked layers
    # are integrated exactly by the sign-chain correlation q=2^{-gap}.
    stop_layers=sorted({1}|{j for j in range(1,len(en)) if en[j]},reverse=True)
    cuts={F(0),F(1)}
    for j in stop_layers:
        for endpoint in en[j]:
            c=abs(endpoint-y)/R[j]
            if 0<c<1:cuts.add(c)
        for endpoint in H.knots:
            c=abs(endpoint-y)/R[j]
            if 0<c<1:cuts.add(c)
    cuts=sorted(cuts);hvalue=mass=absorbed=F(0);stop_law=defaultdict(F);eta_events=defaultdict(F)
    def endpoint_density(j,sign,lo,hi,prob):
        if not prob:return
        l,h=sorted((y+sign*R[j]*lo,y+sign*R[j]*hi));density=prob/R[j]
        eta_events[l]+=density;eta_events[h]-=density
    for lo,hi in zip(cuts,cuts[1:]):
        c=(lo+hi)/2;plus=minus=F(1,2)
        for index,j in enumerate(stop_layers):
            if j==1:
                hvalue+=(hi-lo)*(plus*H(y+R[j]*c)+minus*H(y-R[j]*c))
                endpoint_density(j,1,lo,hi,plus);endpoint_density(j,-1,lo,hi,minus)
                mass+=(hi-lo)*(plus+minus);stop_law[j]+=(hi-lo)*(plus+minus);break
            if inside(y+R[j]*c,en[j]):
                hvalue+=(hi-lo)*plus*H(y+R[j]*c);mass+=(hi-lo)*plus
                endpoint_density(j,1,lo,hi,plus)
                absorbed+=(hi-lo)*plus;stop_law[j]+=(hi-lo)*plus;plus=F(0)
            if inside(y-R[j]*c,en[j]):
                hvalue+=(hi-lo)*minus*H(y-R[j]*c);mass+=(hi-lo)*minus
                endpoint_density(j,-1,lo,hi,minus)
                absorbed+=(hi-lo)*minus;stop_law[j]+=(hi-lo)*minus;minus=F(0)
            next_j=stop_layers[index+1];q=F(1,2**(j-next_j))
            same=(1+q)/2;flip=(1-q)/2
            plus,minus=same*plus+flip*minus,flip*plus+same*minus
    assert mass==1==sum(stop_law.values(),F(0))
    return hvalue,absorbed,len(cuts)-1,stop_law,eta_events


def run_scope(atoms,alpha,scope,ratio):
    saved,en,R,H=observer(atoms,alpha,scope,ratio);etaH=hit=F(0);records=[];stops=defaultdict(F);eta_events=defaultdict(F)
    for y,w in atoms:
        value,prob,cells,law,events=stopped_source(y,en,R,H);etaH+=w*value;hit+=w*prob
        for x,density in events.items():eta_events[x]+=w*density
        for j,v in law.items():stops[j]+=w*v
        records.append(dict(source=str(y),mass=str(w),H_at_source=str(H(y)),conditional_stopped_H=str(value),
          conditional_defect=str(H(y)-value),entrance_absorption_before_final1=str(prob),exact_c_cells=cells))
    a=F(saved['entrance_a']);b=F(saved['terminal_b']);V=F(saved['eligible_volume']);T=F(saved['T_exact'])
    knots=sorted(x for x,v in eta_events.items() if v);eta=[];density=F(0)
    for lo,hi in zip(knots,knots[1:]):
        density+=eta_events[lo];assert 0<=density<=2*a
        if density:eta.append((lo,hi,density))
    assert density+eta_events[knots[-1]]==0
    assert sum(((hi-lo)*d for lo,hi,d in eta),F(0))==1
    assert sum(((hi*hi-lo*lo)*d/2 for lo,hi,d in eta),F(0))==sum((y*w for y,w in atoms),F(0))
    assert sum((d*(H.primitive(hi)-H.primitive(lo)) for lo,hi,d in eta),F(0))==etaH
    free=2*a*F(saved['H_positive_integral_exact']);cap,t,trials=H.cap_upper(2*a)
    coarseH=sum((w*(H.primitive(y+R[1])-H.primitive(y-R[1]))/(2*R[1]) for y,w in atoms),F(0))
    assert coarseH<=free
    assert etaH<=cap<=free<=2*a*V
    if ratio==4:assert etaH<=F(2,3)*T and T-etaH>=T/3
    assert sum(stops.values(),F(0))==1
    saved.update(stopped_H_integral_exact=str(etaH),source_stop_defect_exact=str(T-etaH),
      full_coarse_H_integral_exact=str(coarseH),full_coarse_defect_exact=str(T-coarseH),
      stop_to_coarse_signed_continuation_exact=str(coarseH-etaH),
      exact_continue_decomposition_matches=(T-etaH)==(T-coarseH)+(coarseH-etaH),
      stopped_H_over_T_exact=str(etaH/T) if T else None,defect_over_T_exact=str((T-etaH)/T) if T else None,
      cap_2a_exact=str(2*a),free_Hplus_upper_exact=str(free),selected_cap_upper_exact=str(cap),cap_trial_winning_t=str(t),
      cap_trials=[dict(t=str(t),upper=str(v)) for t,v in trials],cap_trials_are_valid_upper_bounds_not_claimed_optimum=True,
      geometric_free_upper_exact=str(2*a*V),geometric_absorption_coefficient_exact=str(2*a/(b-a)),
      source_ledger=records,complete_stop_mass_exact='1',entrance_absorption_before_final1_exact=str(hit),
      stop_layer_law=[dict(j=j,mass=str(v)) for j,v in sorted(stops.items()) if v],
      complete_eta_density=[dict(lo=str(lo),hi=str(hi),density=str(d)) for lo,hi,d in eta],
      eta_density_peak_exact=str(max(d for lo,hi,d in eta)),eta_density_mass_exact='1',
      eta_first_moment_matches_source=True,eta_density_H_integral_matches_radial_integral=True,
      complete_source_radial_event_cells=sum(v['exact_c_cells'] for v in records),
      displays=dict(T=float(T),etaH=float(etaH),defect=float(T-etaH),etaH_over_T=float(etaH/T) if T else None,
        cap_upper=float(cap),free_Hplus=float(free),coarseH=float(coarseH),signed_continuation=float(coarseH-etaH)))
    return saved


def main():
    start=time.monotonic();root=Path(__file__).resolve().parents[2]
    old=json.loads((root/'output/general_input_20261003/stage20_source_averaged_probe.json').read_text());rounds=[]
    for rd in old['rounds']:
        case=next(c for c in rd['cases'] if c['family']=='four_digital_clusters')
        atoms=sorted((F(v['location']),F(v['mass'])) for v in case['atoms']);scopes=[]
        for ratio in [2,4]:
            for scope in ['full','band']:
                result=run_scope(atoms,F(1),scope,ratio);scopes.append(result)
                print(json.dumps(dict(round=rd['round'],N=len(atoms),scope=scope,b_over_a=ratio,**result['displays'])),flush=True)
        rounds.append(dict(round=rd['round'],source=case['atoms'],parameter=case['parameter'],scopes=scopes))
    out=dict(status='passed',dimension=1,budget='Three frozen digital sources; two entrance parameter choices, each with original full/band, twelve exact scopes.',
      arithmetic='Fraction exact full observer, all sources, common-c Brownian radial/sign stopping and piecewise-linear H integrals; displays only are floats.',
      rounds=rounds,elapsed_seconds=time.monotonic()-start,main_uniform_weak_bound_proved=False,
      source_stop_defect_uniform_bound_proved=False,no_source_or_observer_sampling=True)
    path=root/'output/general_input_20261003/stage24_clock_replacement_probe.json';path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status='passed',seconds=out['elapsed_seconds'])),flush=True)


if __name__=='__main__':main()
