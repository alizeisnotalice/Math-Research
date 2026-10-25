"""Fixed twelve-block parity, signed-radius and raw-endpoint checks. No union search."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json
import sympy as sy
from vg8_parity import split_block,unpack_key,original_to_super,super_to_original,finite_mass,radius_descriptor,compile_rational_receiver,old_residue_keys

P=Path(__file__).resolve().parent
CONFIG_SHA='05b2f7804dfdf8d16749a844fe7ed0361379e1ac02ecfd1f714763f398ffb551'


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    return x


def atoms(lo,hi):
    return list(product(*(range(a,b+1) for a,b in zip(lo,hi))))


def main():
    config_path=P/'vg8_blocks.json'
    assert hashlib.sha256(config_path.read_bytes()).hexdigest()==CONFIG_SHA
    config=json.loads(config_path.read_text());eta=F(config['eta']);epsilons=list(map(F,config['finite_mass_parameters']))
    rows=[];total_atoms=0
    for spec in config['blocks']:
        lo=spec['lo'];hi=spec['hi'];original=atoms(lo,hi);split=split_block(lo,hi)
        kl,kh,counts=unpack_key(split['key'])
        generated={(p,k) for p in (0,1) for k in atoms(kl[p],kh[p])}
        mapped={original_to_super(k) for k in original}
        assert generated==mapped and len(generated)==len(original)
        assert all(super_to_original(p,k) in original for p,k in generated)
        assert counts==(sum(k[0]%2==0 for k in original),sum(k[0]%2==1 for k in original))
        length=hi[0]-lo[0]+1;transverse=len(original)//length
        expected_delta=0 if length%2==0 else (1 if lo[0]%2==0 else -1)*transverse
        assert split['delta']==expected_delta
        residue0,residue1=old_residue_keys(lo,hi)
        assert residue0==split and residue1['N_even']==split['N_odd'] and residue1['N_odd']==split['N_even']
        masses=[]
        for e in epsilons:
            brute=sum((1+e if k[0]%2==0 else 1-e for k in original),F(0))
            assert finite_mass(split['key'],e)==brute==split['N']+e*split['delta']
            assert finite_mass(residue1['key'],e)==finite_mass(residue0['key'],-e)
            assert (1+e)+(1-e)==2
            masses.append(brute)
        translated={original_to_super((k[0]+2,)+k[1:]) for k in original}
        expected={(p,(k[0]+1,)+k[1:]) for p,k in mapped}
        assert translated==expected
        assert all(F(2)*(k[0]+F(p,2))==super_to_original(p,k)[0] for p,k in mapped)
        rows.append(dict(lo=lo,hi=hi,key=split['key'],counts=counts,N=split['N'],delta=split['delta'],masses=masses))
        total_atoms+=len(original)
    assert len(rows)==12
    N,D,t=sy.symbols('N D t',positive=True)
    r=(N+t*D)**sy.Rational(1,4)/2
    assert sy.simplify(sy.diff(r,t).subs(t,0)-D/(8*N**sy.Rational(3,4)))==0
    assert sy.simplify(-sy.diff(r,t,2).subs(t,0)/2-3*D**2/(64*N**sy.Rational(7,4)))==0
    examples=[]
    compiled=[]
    for row in rows[-2:]:
        info=radius_descriptor(row['key'],1,eta)
        compiled_row=compile_rational_receiver(row['key'],list(map(F,config['base_periods_P'])),1,
            list(map(F,config['relative_P_shift_direction_h'])),list(map(F,config['period_direction'])),eta)
        sign=1 if row['delta']>0 else -1
        assert compiled_row['r0']==F(3,2) and compiled_row['r1']==sign*F(1,8)
        expected_y=(F(1,8),F(1,8)) if sign>0 else (F(3,8),F(-1,8))
        assert (compiled_row['lower'][1].velocity,compiled_row['upper'][1].velocity)==expected_y
        examples.append(dict(delta=row['delta'],lambda_r=info['radius_velocity_factor'],
            beta=info['radius_remainder_factor'],r0=compiled_row['r0'],r1=compiled_row['r1'],
            y_velocity=expected_y,complete_mass_slope=info['complete_mass_slope']))
        compiled.append(compiled_row)
    # One finite rational fourth-root check, deliberately on the decreasing-mass block.
    q=F(4095,4096);tau=3*(1-q**4);rtrue=F(3,2)*q;c=compiled[0]
    assert 0<tau<eta and rtrue**4==finite_mass(rows[-2]['key'],tau)/16
    remainder=rtrue-c['r0']-c['r1']*tau
    Br=c['r0']*c['descriptor']['radius_remainder_factor']
    assert remainder<0 and abs(remainder)<=Br*tau**2
    lo,hi,counts=unpack_key(rows[-2]['key']);u=c['super_shift_velocity'];basephase=[(F(0),)*4,(F(1,2),F(0),F(0),F(0))]
    finite_errors=[]
    for j in range(4):
        pos_hi=[hi[p][j]+basephase[p][j]+(tau*u[j] if p else 0) for p in (0,1)]
        pos_lo=[lo[p][j]+basephase[p][j]+(tau*u[j] if p else 0) for p in (0,1)]
        actual=(max(pos_hi)-rtrue/c['Q'][j],min(pos_lo)+rtrue/c['Q'][j])
        for side,value in zip(('lower','upper'),actual):
            e=c[side][j];error=abs(value-e.base-e.velocity*tau)
            assert error<=e.remainder*tau**2
            finite_errors.append(error)
    result=dict(status='fixed parity/signed-radius/endpoint checks PASS; no R4 union recomputation',
        configuration_sha256=CONFIG_SHA,block_count=12,enumerated_atoms=total_atoms,epsilon_values=epsilons,
        rows=rows,N81_examples=examples,
        finite_negative_mass_example=dict(tau=tau,root_ratio=q,rtrue=rtrue,remainder=remainder,
            Br=Br,bound=Br*tau**2,endpoint_absolute_errors=finite_errors),
        source_sha256={name:hashlib.sha256((P/name).read_bytes()).hexdigest() for name in ['vg8_parity.py','vg8_check.py','vg6_bounds.py']})
    (P/'vg8_results.json').write_text(json.dumps(encode(result),indent=2)+'\n')
    print(json.dumps(encode({k:v for k,v in result.items() if k!='rows'}),indent=2))


if __name__=='__main__':main()
