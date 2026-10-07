"""Exact full-observer matched-height ball residuals, no fixed-H substitution.
Old ratio4 frozen scopes plus three new ratio8 rational inputs.
Integration of Phi' endpoint flux uses Phi endpoints, not a height grid.
"""
from fractions import Fraction as F
from pathlib import Path
from bisect import bisect_right
import json,sys,time,hashlib
sys.dont_write_bytecode=True
from stage24_clock_replacement_probe import run_scope
from stage25_stopped_energy_probe import MomentMeasure
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output/general_input_20261003/stage29_matched_height_probe.json'
DELTA=F(1,16)
def phi(v,t):
    assert v>=0 and t>0
    return v-v*v/(2*t) if v<t else t/2
def beta(v,t):return max(F(0),1-v/t)
def measure_mass(P,x):return P.moments[0][bisect_right(P.points,x)]
def audit(atoms,saved,label):
    alpha=F(saved['alpha']);ratio=saved['b_over_a'];a=F(saved['entrance_a']);c=2*a
    P=MomentMeasure(atoms=atoms)
    eta=MomentMeasure(cells=[tuple(F(v[k]) for k in ('lo','hi','density')) for v in saved['complete_eta_density']])
    assert P.total()==eta.total()==1 and P.total(1)==eta.total(1)
    u=lambda x:eta.left_distance(x)-P.left_distance(x)
    du=lambda x:eta.prefix(x)[0]-measure_mass(P,x)
    D=saved['D'];p=[F(1,2**j) for j in range(D+1)];R=[pp/(2*a) for pp in p]
    E=[(F(l),F(h)) for l,h in saved['original_intervals']]
    cuts=sorted({x for l,h in E for x in (l,h)}|{y+sign*r for y,w in atoms for r in R[1:] for sign in (-1,1)})
    endpoints=[x for l,h in E for x in (l,h)]
    inside=lambda x:bisect_right(endpoints,x)%2==1
    rows=[];VP=Veta=Vphi=Vpi=Vres=volume=eligible=J1=F(0);cap=DELTA*alpha+F(3,2)*c
    split={k:dict(volume=F(0),source_average=F(0),matched_Pi_average=F(0),matched_signed_residual=F(0)) for k in ('J1','eligible_J_at_least_2')}
    for l,h in zip(cuts,cuts[1:]):
        x=(l+h)/2
        if not inside(x):continue
        masses=[measure_mass(P,x+r)-measure_mass(P,x-r) for r in R]
        J=next(j for j in range(1,D+1) if masses[j]>p[j])
        K=next(j for j in range(1,D+1) if masses[j]>ratio*p[j])
        assert J<=K
        r=R[K];t=DELTA*alpha*r*r;length=h-l;gP=masses[K]/(2*r)
        assert alpha/2<gP<=alpha
        emass_integral=eta.left_distance(h+r)-eta.left_distance(l+r)-eta.left_distance(h-r)+eta.left_distance(l-r)
        corner_points=(l+r,h+r,l-r,h-r);potentials=tuple(u(z) for z in corner_points)
        assert all(v>=0 for v in potentials)
        phis=tuple(phi(v,t) for v in potentials)
        phi_change=phis[1]-phis[0]-phis[3]+phis[2]
        piavg_integral=(emass_integral-phi_change)/(2*r)
        pavg_integral=length*gP;res=pavg_integral-piavg_integral
        assert 0<=piavg_integral<=cap*length
        assert res>=(alpha/2-cap)*length
        # Strictly inside the original cell, select a rational center whose
        # ball endpoints avoid atomic sources. This is only a point test.
        for fraction in (F(1,2),F(3,7),F(2,5),F(1,3)):
            test=l+(h-l)*fraction
            if all(test+sign*r!=y for y,w in atoms for sign in (-1,1)):break
        else:raise AssertionError('No generic test center')
        ep=eta.prefix(test+r)[0]-eta.prefix(test-r)[0]
        flux=beta(u(test+r),t)*du(test+r)-beta(u(test-r),t)*du(test-r)
        assert u(test+r)>=0 and u(test-r)>=0
        pointpi=(ep-flux)/(2*r)
        assert 0<=pointpi<=cap
        volume+=length
        if J==1:J1+=length
        else:eligible+=length
        VP+=pavg_integral;Veta+=emass_integral/(2*r);Vphi+=phi_change/(2*r);Vpi+=piavg_integral;Vres+=res
        group=split['J1' if J==1 else 'eligible_J_at_least_2']
        for key,value in [('volume',length),('source_average',pavg_integral),('matched_Pi_average',piavg_integral),('matched_signed_residual',res)]:group[key]+=value
        rows.append(dict(lo=str(l),hi=str(h),J=J,K=K,radius=str(r),matched_height=str(t),
          source_ball_mass=str(masses[K]),source_ball_average=str(gP),
          corner_points=list(map(str,corner_points)),corner_u=list(map(str,potentials)),corner_Phi=list(map(str,phis)),
          eta_ball_mass_spatial_integral=str(emass_integral),Phi_endpoint_spatial_flux=str(phi_change),
          source_average_spatial_integral=str(pavg_integral),matched_Pi_average_spatial_integral=str(piavg_integral),
          matched_signed_residual_spatial_integral=str(res),
          rational_point_certificate=dict(center=str(test),eta_ball_mass=str(ep),Phi_derivative_flux=str(flux),
          matched_Pi_ball_average=str(pointpi),point_cap_verified=True),
          exact_whole_cell_integral_cap_verified=True))
    assert volume==F(saved['original_volume']) and eligible==F(saved['eligible_volume']) and J1==F(saved['J1_volume'])
    assert Vpi==Veta-Vphi and Vres==VP-Vpi
    assert Vpi<=cap*volume and Vres>=(alpha/2-cap)*volume
    return dict(label=label,scope=saved['scope'],b_over_a=ratio,alpha=str(alpha),delta=str(DELTA),
      complete_actual_observer_stop=saved,original_full_volume=str(volume),eligible_volume=str(eligible),J1_volume=str(J1),
      complete_original_cells=rows,all_original_observers_covered=True,
      cap_coefficient_over_alpha=str(cap/alpha),condenser_cap=str(cap),
      full_spatial_integrals=dict(original_source_average=str(VP),terminal_eta_average=str(Veta),
        clipped_Phi_endpoint_flux=str(Vphi),matched_Pi_average=str(Vpi),matched_signed_residual=str(Vres)),
      residual_lower_bound=str((alpha/2-cap)*volume),
      separate_J1_and_eligible_integrals={k:{key:str(value) for key,value in row.items()} for k,row in split.items()},
      display=dict(volume=float(volume),source_integral=float(VP),matched_Pi_integral=float(Vpi),
        matched_residual=float(Vres),residual_over_alpha_volume=float(Vres/(alpha*volume))),
      all_matched_cell_integrals_exact=True,point_tests_do_not_certify_every_center_cap=True,
      no_quadratic_root_partition_needed_for_signed_integral=True,
      original_H_not_used_as_a_single_height_pairing=True)
def main():
    start=time.monotonic()
    raw=(ROOT/'output/general_input_20261003/stage24_clock_replacement_probe.json').read_bytes();old=json.loads(raw)
    oldrounds=[]
    for rd in old['rounds']:
        atoms=[(F(v['location']),F(v['mass'])) for v in rd['source']];scopes=[]
        for s in rd['scopes']:
            if s['b_over_a']!=4:continue
            got=audit(atoms,s,f"frozen-N{len(atoms)}-{s['scope']}");scopes.append(got)
            print(got['label'],got['display'],flush=True)
        oldrounds.append(dict(round=rd['round'],source=rd['source'],scopes=scopes))
    cases=[
      ('unequal_near_far',F(1),['0','2/25','2/5','8'],[1,2,3,4]),
      ('geometric_small_cluster',F(3,2),['0','1/64','1/8','1/2','2','6'],[1,2,4,3,2,3]),
      ('fork_and_far',F(5,4),['0','1/128','1/32','1/8','1/2','7/4','3','7'],[1,3,2,5,4,3,2,4])]
    newrounds=[]
    for rd,(name,alpha,xs,ws) in enumerate(cases,1):
        atoms=sorted((F(x),F(w,sum(ws))) for x,w in zip(xs,ws));scopes=[]
        for scope in ('full','band'):
            saved=run_scope(atoms,alpha,scope,8)
            got=audit(atoms,saved,f"ratio8-N{len(atoms)}-{scope}");scopes.append(got)
            print(got['label'],got['display'],flush=True)
        newrounds.append(dict(round=rd,family=name,source=[dict(location=str(x),mass=str(w)) for x,w in atoms],scopes=scopes))
    out=dict(status='passed',dimension=1,delta=str(DELTA),
      method='Exact Fraction endpoint Phi identity integrates every actual original observer cell; rational generic-center ball certificates separately saved.',
      frozen_ratio4_input=dict(file='stage24_clock_replacement_probe.json',sha256=hashlib.sha256(raw).hexdigest()),
      frozen_ratio4_rounds=oldrounds,new_ratio8_rounds=newrounds,
      counts=dict(old_full_band_scopes=6,new_full_band_scopes=6,new_inputs=3,
        complete_original_observer_cells=sum(len(s['complete_original_cells']) for r in oldrounds+newrounds for s in r['scopes'])),
      no_source_observer_or_height_sampling_for_full_spatial_integral=True,
      main_weak_theorem_proved=False,universal_volume_budget_proved_by_experiments=False,
      elapsed_seconds=time.monotonic()-start)
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('completed',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
