"""Independent Fraction audit via power conversion and affine substitution.

No de Casteljau/repeated degree-elevation implementation is used here.
No Arb, integral, function point, AS/Psi model or LP is evaluated.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb
import json
from pathlib import Path
from time import monotonic

ROOT = Path(__file__).resolve().parent
SOURCE_SHA = 'fd16ccf0309a4ae4627049ff730b2bb233be3948013d25cd514824febcc107b2'
MODEL_SHA = '8b614d16e19e8dddd9cb16793595b99aff3acf680795710c50ced9c3f76c108e'
COVER_SHA = 'e78b2bf4e6f5e00974a180cb85dfa9044ed42939930c4add48078e6c4064c8e7'
BASE_SHA = '21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
PHYSICAL = ('6677/4096','3339/2048','161/4096','81/2048')
LOCAL = ('21/128','11/64','33/128','17/64')


def bernstein_to_power(B):
    """Expand C(3,i)t^i(1-t)^(3-i) in each coordinate."""
    P = [[Q(0) for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            for k in range(i,4):
                for l in range(j,4):
                    P[k][l] += B[i][j]*comb(3,i)*comb(3-i,k-i)*(-1)**(k-i)*comb(3,j)*comb(3-j,l-j)*(-1)**(l-j)
    return P


def affine_substitution(P, a, b, c, d):
    """Substitute old coordinates t=a+(b-a)x, s=c+(d-c)y."""
    out = [[Q(0) for _ in range(4)] for _ in range(4)]
    for k in range(4):
        for l in range(4):
            for i in range(k+1):
                for j in range(l+1):
                    out[i][j] += P[k][l]*comb(k,i)*a**(k-i)*(b-a)**i*comb(l,j)*c**(l-j)*(d-c)**j
    return out


def power_to_bernstein(P):
    """Use t^k=sum_{i>=k} C(i,k)/C(3,k) B_i^3(t)."""
    return [[sum((P[k][l]*Q(comb(i,k),comb(3,k))*Q(comb(j,l),comb(3,l))
        for k in range(i+1) for l in range(j+1)),Q(0)) for j in range(4)] for i in range(4)]


def restrict(B, domain):
    return power_to_bernstein(affine_substitution(bernstein_to_power(B),*domain))


def matrix(rows):
    assert len(rows)==4 and all(len(row)==4 for row in rows)
    return [list(map(Q,row)) for row in rows]


def main():
    output = ROOT/'saved_verification.json'
    if output.exists():
        raise SystemExit('Existing independent audit: no repeat or overwrite.')
    start = monotonic()
    def guard():
        if monotonic()-start > 90:
            raise TimeoutError('90 second exact saved-value audit limit')
    def read(path):
        guard();return json.loads(path.read_text())
    def digest(path):
        guard();return sha256(path.read_bytes()).hexdigest()
    result, proof, run = (read(ROOT/name) for name in ('result.json','quadrants.json','quadrants.run.json'))
    assert run['status']=='terminal' and run['exit_code']==0 and 0<=run['elapsed_seconds']<=60
    assert run['source_sha256_before']==run['source_sha256_after']==result['source_sha256']==digest(ROOT/'quadrants_once.py')==SOURCE_SHA
    assert run['log_sha256']==digest(ROOT/'quadrants.log') and result['quadrants_sha256']==digest(ROOT/'quadrants.json')
    expected = {'round22/models.local/node_00663.json':MODEL_SHA,
        'round22/leaf_cover.json':COVER_SHA,'round01/scaled_certificate.json':BASE_SHA}
    assert len(result['inputs'])==3 and {v['file']:v['sha256'] for v in result['inputs']}==expected
    for relative,wanted in expected.items(): assert digest(ROOT.parent/relative)==wanted
    old = read(ROOT.parent/'round22/models.local/node_00663.json')
    cover = read(ROOT.parent/'round22/leaf_cover.json');base = read(ROOT.parent/'round01/scaled_certificate.json')
    node = cover['nodes'][663]
    assert node['status']=='unresolved_depth4' and node['depth']==4 and node['model_sha256']==MODEL_SHA
    assert node['parent']==old['parent']==[52,1] and node['id']==old['node_id']==proof['node_id']==663
    assert tuple(node['physical'])==tuple(old['physical'])==tuple(proof['physical'])==PHYSICAL
    assert tuple(node['local'])==tuple(old['local'])==tuple(proof['local'])==LOCAL
    assert proof['source_model_sha256']==MODEL_SHA and proof['source_cover_sha256']==COVER_SHA and proof['baseline_sha256']==BASE_SHA
    assert proof['source_tree_depth']==result['source_tree_depth']==4
    assert proof['source_tree_changed'] is False and result['source_tree_changed'] is False
    assert proof['proof_partition']=='four_center_sign_quadrants'
    assert (Q(base['alpha']),Q(base['R']),Q(base['h']),base['N'])==(Q(47,40),Q(33,20),Q(1,32),80)
    c1,h,R = Q(base['c']),Q(base['h']),Q(base['R']);assert c1>0
    u0,u1,v0,v1 = map(Q,PHYSICAL);x0,x1,y0,y1 = map(Q,LOCAL)
    uc,vc = (u0+u1)/2,(v0+v1)/2;xc,yc = (x0+x1)/2,(y0+y1)/2
    du,dv = (u1-u0)/2,(v1-v0)/2
    assert u0>0 and v0>0 and u0+v0>R
    assert (u0,u1,v0,v1)==(h*(52+x0),h*(52+x1),h*(1+y0),h*(1+y1))
    assert list(map(Q,proof['physical_center']))==[uc,vc] and list(map(Q,proof['local_center']))==[xc,yc]
    psi = old['Psi'];es = Q(old['AS']['remainder_upper']);muv = Q(psi['Muv_upper']);D = Q(psi['one_sided_lower_deduction'])
    assert psi['one_sided_only'] and psi['unscaled_by_c'] and psi['Bernstein_degree']==[3,3]
    assert es>=0 and muv>=0 and D==muv*du*dv and Q(psi['du'])==du and Q(psi['dv'])==dv
    assert Q(proof['AS_remainder_upper'])==es and Q(proof['Muv_upper'])==muv and Q(proof['old_uniform_D'])==D
    absolute_uppers = []
    assert len(psi['mixed_derivatives'])==2
    for br,(a,b) in zip(psi['mixed_derivatives'],((u0,v0),(v0,u0))):
        lo,hi,au = map(Q,(br['lower'],br['upper'],br['absolute_upper']))
        assert (Q(br['a']),Q(br['b']),br['kind'])==(a,b,'Hab')
        assert lo<=hi and au>=max(abs(lo),abs(hi));absolute_uppers.append(au)
    assert sum(absolute_uppers,Q(0))==muv
    arrays = {name:[[Q(0) for _ in range(4)] for _ in range(4)] for name in ('T1','AS','old_c','old_s','old_slack')}
    assert len(psi['bernstein'])==len(old['AS']['bernstein'])==len(old['rows'])==16
    for k in range(16):
        i,j = divmod(k,4);pb,sb,row = psi['bernstein'][k],old['AS']['bernstein'][k],old['rows'][k]
        assert pb['index']==sb['index']==row['bernstein_index']==[i,j]
        pl,pu,sl,su = map(Q,(pb['lower'],pb['upper'],sb['lower'],sb['upper']))
        assert pl<=pu and sl<=su
        assert (Q(row['c_column']),Q(row['s_column']),Q(row['baseline_slack']))==(pl-D,sl-es,c1*(pl-D)+sl-es)
        for name,value in (('T1',pl),('AS',sl),('old_c',pl-D),('old_s',sl-es),('old_slack',c1*(pl-D)+sl-es)):
            arrays[name][i][j]=value
    assert matrix(proof['source_T1_lower'])==arrays['T1'] and matrix(proof['source_AS_lower'])==arrays['AS']
    # Independent algorithm self-checks only exact polynomial identities.
    for B in arrays.values(): assert power_to_bernstein(bernstein_to_power(B))==B
    assert len(proof['pieces'])==4
    physical_rects=[];all_rows=[];all_gains=[];piece_acceptance=[]
    halves = ((Q(0),Q(1,2)),(Q(1,2),Q(1)))
    for pid,(iu,iv) in enumerate(((0,0),(1,0),(0,1),(1,1))):
        guard();piece = proof['pieces'][pid]
        assert piece['piece_id']==pid and piece['proof_piece'] is True and piece['source_tree_depth']==4
        assert (piece['u_sign'],piece['v_sign'])==('-' if iu==0 else '+','-' if iv==0 else '+')
        assert piece['same_sign']==(iu==iv)
        xa,xb = halves[iu];ya,yb = halves[iv];domain = (xa,xb,ya,yb)
        a,b = u0+(u1-u0)*xa,u0+(u1-u0)*xb
        c,d = v0+(v1-v0)*ya,v0+(v1-v0)*yb
        assert tuple(map(Q,piece['physical']))==(a,b,c,d)
        assert tuple(map(Q,piece['local']))==(x0+(x1-x0)*xa,x0+(x1-x0)*xb,y0+(y1-y0)*ya,y0+(y1-y0)*yb)
        physical_rects.append((a,b,c,d))
        restricted = {name:restrict(B,domain) for name,B in arrays.items()}
        assert matrix(piece['restricted_T1_lower'])==restricted['T1'] and matrix(piece['restricted_AS_lower'])==restricted['AS']
        # Build correction directly in the new piece POWER basis, not by elevation.
        correction_power = [[Q(0) for _ in range(4)] for _ in range(4)]
        if iu!=iv:
            correction_power[0][0]=muv*(a-uc)*(c-vc)
            correction_power[1][0]=muv*(b-a)*(c-vc)
            correction_power[0][1]=muv*(a-uc)*(d-c)
            correction_power[1][1]=muv*(b-a)*(d-c)
        correction = power_to_bernstein(correction_power)
        assert matrix(piece['bilinear_Bernstein_degree3'])==correction
        assert len(piece['bilinear_Bernstein_degree1'])==2
        corners = [[sum((correction_power[k][l]*Q(i)**k*Q(j)**l for k in range(4) for l in range(4)),Q(0)) for j in range(2)] for i in range(2)]
        assert [list(map(Q,line)) for line in piece['bilinear_Bernstein_degree1']]==corners
        assert len(piece['rows'])==16
        local_slacks=[]
        for k,row in enumerate(piece['rows']):
            guard();i,j=divmod(k,4);assert row['bernstein_index']==[i,j]
            q=correction[i][j];cc=restricted['T1'][i][j]+q;ss=restricted['AS'][i][j]-es
            rho=c1*cc+ss;oc,os,old_rho=(restricted[name][i][j] for name in ('old_c','old_s','old_slack'))
            assert oc==restricted['T1'][i][j]-D and os==ss and old_rho==c1*oc+os
            gain=rho-old_rho;assert gain==c1*(D+q) and gain>=0
            assert tuple(Q(row[name]) for name in ('c_column','s_column','baseline_slack','restricted_old_c_column','restricted_old_s_column','restricted_old_slack','gain'))==(cc,ss,rho,oc,os,old_rho,gain)
            all_rows.append((rho,pid,i,j));all_gains.append(gain);local_slacks.append(rho)
        accepted=all(v>0 for v in local_slacks)
        assert piece['accepted']==accepted;piece_acceptance.append(accepted)
    # Independent physical slab coverage: four pieces partition the ORIGINAL node.
    xs = sorted({u0,u1}|{x for rect in physical_rects for x in rect[:2]})
    for left,right in zip(xs,xs[1:]):
        cursor=v0
        for bottom,top in sorted((c,d) for a,b,c,d in physical_rects if a<=left and right<=b):
            assert bottom==cursor;cursor=top
        assert cursor==v1
    assert sum(((b-a)*(d-c) for a,b,c,d in physical_rects),Q(0))==(u1-u0)*(v1-v0)
    minimum=min(all_rows)
    assert result['status']=='completed' and result['accepted']==all(piece_acceptance)
    assert result['pieces_started']==result['pieces_completed']==4 and result['rows_started']==result['rows_completed']==64
    assert result['in_flight'] is None and Q(result['minimum_slack'])==minimum[0]
    assert result['minimum_piece']==minimum[1] and result['minimum_index']==[minimum[2],minimum[3]]
    assert Q(result['minimum_gain'])==min(all_gains)
    assert all(result[name]==0 for name in ('new_Arb_calls','new_integrals','new_AS','new_Psi_models','new_points','LP_calls'))
    assert digest(ROOT.parent/'round22/leaf_cover.json')==COVER_SHA
    answer=dict(status='passed_independent_power_substitution_audit',
        algorithm='Bernstein to power, exact affine substitution, power to Bernstein; no de Casteljau calls.',
        scope='Saved lower polynomials, 4 fixed proof pieces/64 rows and original-node partition only; no new function/model/integral evaluation.',
        pieces_checked=4,rows_checked=64,physical_slabs_checked=len(xs)-1,accepted=all(piece_acceptance),
        minimum_slack=str(minimum[0]),minimum_piece=minimum[1],minimum_index=[minimum[2],minimum[3]],
        minimum_gain=str(min(all_gains)),strictly_positive_gains=sum(v>0 for v in all_gains),zero_gains=sum(v==0 for v in all_gains),
        source_tree_unchanged=True,source_tree_depth=4,new_integrals=0,new_Arb_evaluations=0,
        elapsed_seconds=monotonic()-start,result_sha256=digest(ROOT/'result.json'),
        quadrants_sha256=digest(ROOT/'quadrants.json'),run_sha256=digest(ROOT/'quadrants.run.json'),
        verifier_sha256=digest(Path(__file__)))
    out=ROOT/'saved_verification.json';tmp=out.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(answer,indent=2)+'\n');tmp.replace(out)
    print(json.dumps(answer))


if __name__=='__main__':
    main()
