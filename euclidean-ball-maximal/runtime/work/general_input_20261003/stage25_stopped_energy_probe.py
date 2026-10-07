#!/usr/bin/env python3
"""Exact stopped occupation energy and Ito square identity, frozen stage24.
No import of the verifier with top-level execution. Prefix moments integrate
the complete source/stopped measures and compact piecewise-linear H exactly.
"""
from fractions import Fraction as F
from bisect import bisect_right
from pathlib import Path
from decimal import Decimal,localcontext
import json,time,hashlib


class MomentMeasure:
    def __init__(self,atoms=None,cells=None):
        self.atoms=atoms or [];self.cells=cells or []
        self.points=[x for x,w in self.atoms] if self.atoms else [l for l,h,d in self.cells]
        self.moments=[[F(0)] for _ in range(3)]
        for row in self.atoms or self.cells:
            for k in range(3):
                value=row[1]*row[0]**k if self.atoms else row[2]*(row[1]**(k+1)-row[0]**(k+1))/(k+1)
                self.moments[k].append(self.moments[k][-1]+value)
    def prefix(self,x):
        i=bisect_right(self.points,x)
        if self.atoms or i==0:return tuple(v[i] for v in self.moments)
        l,h,d=self.cells[i-1];end=min(x,h)
        return tuple(v[i-1]+d*(end**(k+1)-l**(k+1))/(k+1) for k,v in enumerate(self.moments))
    def total(self,k=0):return self.moments[k][-1]
    def left_distance(self,x):
        m,s,t=self.prefix(x);return x*m-s
    def distance_primitive(self,x):
        m,s,t=self.prefix(x);return (x*x*m-2*x*s+t)/2


class LinearFunction:
    def __init__(self,records):
        self.segments=[(F(v['lo']),F(v['hi']),F(v['Hlo']),F(v['Hhi']),F(v['slope'])) for v in records]
        self.points=[lo for lo,hi,v,w,s in self.segments]
        self.square_prefix=[F(0)];self.jumps=[];previous=F(0)
        for index,(lo,hi,v,w,s) in enumerate(self.segments):
            assert w==v+s*(hi-lo)
            if index:assert lo==self.segments[index-1][1] and v==self.segments[index-1][3]
            d=hi-lo;self.square_prefix.append(self.square_prefix[-1]+v*v*d+v*s*d*d+s*s*d**3/3)
            self.jumps.append((lo,v,s-previous));previous=s
        self.jumps.append((hi,w,-previous));assert self.segments[0][2]==w==0
    def value(self,x):
        i=bisect_right(self.points,x)-1
        if i<0:return F(0)
        lo,hi,v,w,s=self.segments[i]
        return v+s*(x-lo) if x<=hi else F(0)
    def square_primitive(self,x):
        i=bisect_right(self.points,x)-1
        if i<0:return F(0)
        lo,hi,v,w,s=self.segments[i]
        if x>=hi:return self.square_prefix[i+1]
        d=x-lo
        return self.square_prefix[i]+v*v*d+v*s*d*d+s*s*d**3/3


def check_scope(atoms,saved):
    P=MomentMeasure(atoms=atoms)
    eta=MomentMeasure(cells=[(F(v['lo']),F(v['hi']),F(v['density'])) for v in saved['complete_eta_density']])
    assert P.total()==eta.total()==1 and P.total(1)==eta.total(1)
    H=LinearFunction(saved['H_linear_segments'])
    potential=lambda x:eta.left_distance(x)-P.left_distance(x)
    primitive=lambda x:eta.distance_primitive(x)-P.distance_primitive(x)
    Q=F(0);segment_records=[]
    for lo,hi,v,w,slope in H.segments:
        mass=primitive(hi)-primitive(lo);assert mass>=0
        energy=slope*slope*mass;Q+=energy
        segment_records.append(dict(lo=str(lo),hi=str(hi),slope=str(slope),occupation_mass=str(mass),energy=str(energy)))
    S=positive=negative=defect=F(0);knots=[]
    for z,value,jump in H.jumps:
        u=potential(z);assert u>=0
        contribution=u*value*jump;S+=contribution
        positive+=max(contribution,F(0));negative+=max(-contribution,F(0));defect-=u*jump
        knots.append(dict(x=str(z),H=str(value),slope_jump=str(jump),u=str(u),drift_square_contribution=str(contribution)))
    assert S==positive-negative and defect==F(saved['source_stop_defect_exact'])
    EP=sum((w*H.value(y)**2 for y,w in atoms),F(0))
    Eeta=sum((d*(H.square_primitive(hi)-H.square_primitive(lo)) for lo,hi,d in eta.cells),F(0))
    change=Eeta-EP
    assert change==2*Q+2*S
    occupation=(eta.total(2)-P.total(2))/2
    assert occupation>=0
    # Independently split at ALL measure and H breakpoints. u is quadratic
    # between them; integrate u and u^2 with exact polynomial formulas.
    cutset=set(H.points)|{H.segments[-1][1]}|set(P.points)|set(eta.points)|{h for l,h,d in eta.cells}
    left,right=H.segments[0][0],H.segments[-1][1]
    cuts=sorted(x for x in cutset if left<=x<=right);A=B=Qcheck=F(0);quadratic_rows=[]
    for lo,hi in zip(cuts,cuts[1:]):
        mid=(lo+hi)/2;i=bisect_right(H.points,mid)-1;slope=H.segments[i][4]
        if not slope:continue
        k=bisect_right(eta.points,mid)-1
        density=eta.cells[k][2] if k>=0 and mid<eta.cells[k][1] else F(0)
        r0=potential(lo);r1=eta.prefix(lo)[0]-P.prefix(lo)[0];r2=density/2;L=hi-lo
        imass=r0*L+r1*L*L/2+r2*L**3/3
        isquare=r0*r0*L+r0*r1*L*L+(r1*r1+2*r0*r2)*L**3/3+r1*r2*L**4/2+r2*r2*L**5/5
        assert imass>=0 and isquare>=0 and potential(hi)==r0+r1*L+r2*L*L
        A+=slope*slope*isquare;B+=slope*slope*L;Qcheck+=slope*slope*imass
        quadratic_rows.append(dict(lo=str(lo),hi=str(hi),u_coefficients_from_lo=list(map(str,(r0,r1,r2))),
          H_slope=str(slope),integral_u=str(imass),integral_u_squared=str(isquare)))
    assert Qcheck==Q and Q*Q<=A*B
    def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)
    with localcontext() as ctx:
        ctx.prec=80
        optimum=2*dec(Q)+2*dec(A*B).sqrt();best_delta=dec(A/B).sqrt()
    D=defect;T=F(saved['T_exact'])
    return dict(scope=saved['scope'],b_over_a=saved['b_over_a'],D_cutoff=saved['D'],
      original_volume=saved['original_volume'],alphaV=saved['alphaV'],eligible_volume=saved['eligible_volume'],
      entrance_a=saved['entrance_a'],terminal_b=saved['terminal_b'],T_exact=str(T),defect_exact=str(D),
      Q_exact=str(Q),Q_over_defect_squared_exact=str(Q/(D*D)),Q_over_T_squared_exact=str(Q/(T*T)),
      A_u_squared_gradient_energy_exact=str(A),B_unweighted_gradient_energy_exact=str(B),
      weighted_Cauchy_Q_squared_le_AB_exactly_verified=True,
      Q_delta_exact_formula='A/delta+2Q+delta B',optimal_delta_exact_expression='sqrt(A/B)',
      optimal_Q_delta_exact_expression='2Q+2sqrt(AB)',optimal_Q_delta_decimal=str(optimum),
      optimal_delta_decimal=str(best_delta),optimal_Q_delta_lower_exact=str(4*Q),
      decimal_square_roots_are_diagnostics_not_interval_certificates=True,
      source_H_square_expectation_exact=str(EP),stopped_H_square_expectation_exact=str(Eeta),
      square_expectation_difference_exact=str(change),signed_Ito_drift_pairing_exact=str(S),
      drift_pairing_positive_exact=str(positive),drift_pairing_negative_exact=str(negative),
      martingale_bracket_expectation_exact=str(2*Q),Ito_square_identity_exactly_verified=True,
      signed_knot_defect_identity_exactly_verified=True,occupation_mass_exact=str(occupation),
      full_source_and_eta_masses_exact='1',moments_zero_one_match=True,
      energy_segments=segment_records,drift_knots=knots,
      quadratic_occupation_cells=quadratic_rows,
      display=dict(Q=float(Q),defect=float(D),Q_over_defect_squared=float(Q/(D*D)),
        EP_H2=float(EP),Eeta_H2=float(Eeta),square_change=float(change),drift_pairing=float(S),occupation=float(occupation),
        A=float(A),B=float(B),optimal_Q_delta=float(optimum)))


def main():
    start=time.monotonic();root=Path(__file__).resolve().parents[2]
    path=root/'output/general_input_20261003/stage24_clock_replacement_probe.json';raw=path.read_bytes();source=json.loads(raw)
    rounds=[]
    for row in source['rounds']:
        atoms=sorted((F(v['location']),F(v['mass'])) for v in row['source']);scopes=[]
        for saved in row['scopes']:
            result=check_scope(atoms,saved);scopes.append(result)
            print(json.dumps(dict(round=row['round'],N=len(atoms),scope=result['scope'],b_over_a=result['b_over_a'],**result['display'])),flush=True)
        rounds.append(dict(round=row['round'],source=row['source'],parameter=row['parameter'],scopes=scopes))
    out=dict(status='passed',input_file=str(path.relative_to(root)),input_sha256=hashlib.sha256(raw).hexdigest(),
      budget='Three frozen stage24 inputs, twelve complete original full/band old/low-entrance scopes; no new cases or sampling.',
      arithmetic='Fraction prefix moments, exact quadratic occupation primitives, H square primitives and signed slope-jump pairing.',
      dimension=1,rounds=rounds,elapsed_seconds=time.monotonic()-start,
      Ito_square_identity='Eeta H^2−EP H^2=2Q+2 integral u H Hsecond',
      general_uniform_energy_upper_bound_proved=False,main_uniform_weak_bound_proved=False)
    outpath=root/'output/general_input_20261003/stage25_stopped_energy_probe.json';outpath.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status='passed',scopes=12,seconds=out['elapsed_seconds'])),flush=True)


if __name__=='__main__':main()
