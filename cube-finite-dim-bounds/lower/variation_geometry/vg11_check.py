"""Locked <=24 block fixtures; exact four-coset checks, no optimization."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from math import prod
import hashlib,json
import sympy as sy
from vg11_chess import *

P=Path(__file__).resolve().parent
CONFIG_SHA='a9490b77633e191264569fa648970e84585ba5692053b584c64bcbb5e9271baf'


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x


def atoms(lo,hi):return list(product(*(range(a,b+1) for a,b in zip(lo,hi))))


def main():
    path=P/'vg11_blocks.json';assert hashlib.sha256(path.read_bytes()).hexdigest()==CONFIG_SHA
    config=json.loads(path.read_text());eta=F(config['eta']);epsilons=list(map(F,config['epsilon_values']))
    period=list(map(F,config['base_periods_P']));h=list(map(F,config['relative_P_shift_h_direction']))
    kappa=list(map(F,config['period_direction']));alpha=F(config['mass_direction_alpha'])
    rows=[];total=0
    for spec in config['blocks']:
        lo=spec['lo'];hi=spec['hi'];old=atoms(lo,hi);split=split_block(lo,hi)
        kl,kh,counts=unpack_key(split['key'])
        generated={(i,k) for i in range(4) for k in atoms(kl[i],kh[i])}
        mapped={original_to_super(k) for k in old}
        assert generated==mapped and len(generated)==len(old)
        assert all(super_to_original(i,k) in old for i,k in generated)
        brute_counts=tuple(sum(original_to_super(k)[0]==i for k in old) for i in range(4))
        assert counts==brute_counts
        expected_delta=0
        if old:
            lengths=[b-a+1 for a,b in zip(lo,hi)]
            d1=0 if lengths[0]%2==0 else (-1 if lo[0]%2 else 1)
            d2=0 if lengths[1]%2==0 else (-1 if lo[1]%2 else 1)
            expected_delta=d1*d2*prod(lengths[2:])
        assert split['delta']==expected_delta
        masses=[]
        for e in epsilons:
            direct=sum((1+(1 if (k[0]+k[1])%2==0 else -1)*e for k in old),F(0))
            assert finite_mass(split['key'],e)==direct==split['N']+e*split['delta']
            assert sum(1+s*e for s in SIGNS)==4
            masses.append(direct)
        residues=old_residue_keys(lo,hi)
        for (a,b),residue in zip(COSETS,residues):
            sign=(-1)**(a+b)
            assert residue['N']==split['N'] and residue['delta']==sign*split['delta']
            for e in epsilons:assert finite_mass(residue['key'],e)==finite_mass(split['key'],sign*e)
        for axis in (0,1):
            shifted=set()
            for k in old:
                z=list(k);z[axis]+=2;shifted.add(original_to_super(z))
            expected={(i,tuple(v+(1 if j==axis else 0) for j,v in enumerate(k))) for i,k in mapped}
            assert shifted==expected
        for i,k in mapped:
            p,r=COSETS[i];source=super_to_original(i,k)
            assert 2*(F(k[0])+F(p,2))==source[0] and 2*(F(k[1])+F(r,2))==source[1]
        rows.append(dict(lo=lo,hi=hi,key=split['key'],counts=counts,N=split['N'],delta=split['delta'],masses=masses))
        total+=len(old)
    assert len(rows)==24
    N=sy.symbols('N',positive=True);D=sy.symbols('D',real=True);t=sy.symbols('t',real=True)
    rad=(N+D*t)**sy.Rational(1,4)/2
    assert sy.simplify(sy.diff(rad,t).subs(t,0)-D/(8*N**sy.Rational(3,4)))==0
    assert sy.simplify(-sy.diff(rad,t,2).subs(t,0)/2-3*D**2/(64*N**sy.Rational(7,4)))==0
    examples=[];compiled=[]
    for row in rows[12:16]:
        c=compile_rational_receiver(row['key'],period,alpha,h,kappa,eta)
        anchor=compile_rational_receiver(row['key'],period,alpha,h,kappa,eta,anchor=True)
        sign=1 if row['delta']>0 else -1
        assert c['r0']==F(3,2) and c['r1']==sign*F(1,24)
        expected_z=(F(1,24),F(-1,24)) if sign>0 else (F(1,8),F(-1,8))
        assert (c['lower'][2].velocity,c['upper'][2].velocity)==expected_z
        for j in range(4):
            for side in ('lower','upper'):
                assert anchor[side][j].base==c[side][j].base
                assert anchor[side][j].velocity==c[side][j].velocity-c['super_shift_velocity'][j]
                assert anchor[side][j].remainder==c[side][j].remainder
        examples.append(dict(counts=row['counts'],delta=row['delta'],r0=c['r0'],r1=c['r1'],
                             lambda_r=c['descriptor']['radius_velocity_factor'],beta=c['descriptor']['radius_remainder_factor'],z_velocity=expected_z))
        compiled.append(c)
    selected_index=config['finite_geometry_block_index'];selected=rows[selected_index]
    c=compiled[selected_index-12];root=F(config['finite_root_ratio']);tau=9*(1-root**4)
    assert selected['delta']==-9 and 0<tau<eta
    rtrue=F(3,2)*root
    assert rtrue**4==finite_mass(selected['key'],tau)/16
    Br=c['r0']*c['descriptor']['radius_remainder_factor']
    radial_error=rtrue-c['r0']-c['r1']*tau
    assert radial_error<0 and abs(radial_error)<=Br*tau**2
    phases=coset_phases([tau*x for x in h]);u=c['super_shift_velocity']
    kl,kh,counts=unpack_key(selected['key']);lower=[];upper=[];errors=[]
    for j in range(4):
        a=max(kh[i][j]+phases[i][j] for i in range(4))-rtrue/c['Q'][j]
        b=min(kl[i][j]+phases[i][j] for i in range(4))+rtrue/c['Q'][j]
        lower.append(a);upper.append(b);assert b>a
        for side,value in (('lower',a),('upper',b)):
            endpoint=c[side][j];error=abs(value-endpoint.base-endpoint.velocity*tau)
            assert error<=endpoint.remainder*tau**2;errors.append(error)
    volume=prod(b-a for a,b in zip(lower,upper))
    # Finite phase wrapping: phase_i+t_i, atom index_i-t_i must agree exactly.
    additions=[tuple(-(x.numerator//x.denominator) for x in row) for row in phases]
    wrapped=[tuple(x+y for x,y in zip(row,shift)) for row,shift in zip(phases,additions)]
    wrapped_key=rephase_key(selected['key'],additions);wl,wh,wc=unpack_key(wrapped_key)
    assert wc==counts
    for i in range(4):
        for j in range(4):
            assert wl[i][j]+wrapped[i][j]==kl[i][j]+phases[i][j]
            assert wh[i][j]+wrapped[i][j]==kh[i][j]+phases[i][j]
    result=dict(status='24 locked block fixtures PASS; four-coset interface only, no R4 union/search',
        configuration_sha256=CONFIG_SHA,block_count=24,block_atom_occurrences=total,epsilon_values=epsilons,
        rows=rows,N81_examples=examples,
        finite_geometry=dict(block_index=selected_index,tau=tau,rtrue=rtrue,radial_remainder=radial_error,
            Br=Br,radial_bound=Br*tau**2,lower=lower,upper=upper,normalized_receiver_volume=volume,
            complete_mass=4,Q_jacobian=4,endpoint_absolute_errors=errors,wrapped_phase_additions=additions),
        source_sha256={name:hashlib.sha256((P/name).read_bytes()).hexdigest() for name in ['vg11_chess.py','vg11_check.py','vg6_bounds.py']})
    (P/'vg11_results.json').write_text(json.dumps(encode(result),indent=2)+'\n')
    print(json.dumps(encode({k:v for k,v in result.items() if k!='rows'}),indent=2))


if __name__=='__main__':main()
