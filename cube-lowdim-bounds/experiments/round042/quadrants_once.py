#!/usr/bin/env python3
"""Four fixed algebraic proof pieces / 64 exact rows / 60s; no Arb or LP."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parent
SECONDS=60.
MODEL_SHA='8b614d16e19e8dddd9cb16793595b99aff3acf680795710c50ced9c3f76c108e'
COVER_SHA='e78b2bf4e6f5e00974a180cb85dfa9044ed42939930c4add48078e6c4064c8e7'
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
PHYSICAL=['6677/4096','3339/2048','161/4096','81/2048']
LOCAL=['21/128','11/64','33/128','17/64']


def dump(path,value):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n');tmp.replace(path)


def split_half(vector):
    """de Casteljau at exactly 1/2; left/right Bernstein arrays, same degree."""
    work=list(vector);left=[work[0]];right=[work[-1]]
    while len(work)>1:
        work=[(work[k]+work[k+1])/2 for k in range(len(work)-1)]
        left.append(work[0]);right.append(work[-1])
    return left,list(reversed(right))


def tensor_quadrants(matrix):
    """First split u columns, then v rows; keys (u_side,v_side)."""
    ul=[[Q(0) for j in range(4)] for i in range(4)]
    ur=[[Q(0) for j in range(4)] for i in range(4)]
    for j in range(4):
        left,right=split_half([matrix[i][j] for i in range(4)])
        for i in range(4):ul[i][j]=left[i];ur[i][j]=right[i]
    out={}
    for iu,u_matrix in enumerate((ul,ur)):
        vl=[];vr=[]
        for row in u_matrix:
            left,right=split_half(row);vl.append(left);vr.append(right)
        out[(iu,0)]=vl;out[(iu,1)]=vr
    return out


def elevate(vector,degree):
    current=list(vector);n=len(current)-1
    while n<degree:
        new=[current[0]]
        for i in range(1,n+1):new.append(Q(i,n+1)*current[i-1]+Q(n+1-i,n+1)*current[i])
        new.append(current[-1]);current=new;n+=1
    return current


def elevate_bilinear(matrix):
    rows=[elevate(row,3) for row in matrix]
    columns=[elevate([rows[i][j] for i in range(2)],3) for j in range(4)]
    return [[columns[j][i] for j in range(4)] for i in range(4)]


def strings(matrix):return [[str(x) for x in row] for row in matrix]


def main():
    if any((ROOT/name).exists() for name in ('result.json','quadrants.json')):
        raise SystemExit('Existing evidence: no repeated algebraic proof experiment.')
    start=monotonic();result_path=ROOT/'result.json'
    result={'status':'started','limits':{'seconds_hard':60,'pieces':4,'new_rows':64},
        'pieces_started':0,'pieces_completed':0,'rows_started':0,'rows_completed':0,
        'new_Arb_calls':0,'new_integrals':0,'new_AS':0,'new_Psi_models':0,'new_points':0,'LP_calls':0,
        'source_tree_changed':False,'source_tree_depth':4,'in_flight':None}
    def guard():
        if monotonic()-start>=SECONDS:raise TimeoutError('internal_60s_deadline')
    def save():result['elapsed_seconds']=monotonic()-start;dump(result_path,result)
    def verify(relative,wanted):
        guard();path=ROOT.parent/relative;actual=sha256(path.read_bytes()).hexdigest()
        if actual!=wanted:raise ValueError('Frozen input SHA mismatch: '+relative)
        result.setdefault('inputs',[]).append({'file':relative,'sha256':actual});return path
    save();proof=None
    try:
        model_path=verify('round22/models.local/node_00663.json',MODEL_SHA)
        cover_path=verify('round22/leaf_cover.json',COVER_SHA)
        base_path=verify('round01/scaled_certificate.json',BASELINE_SHA)
        old=json.loads(model_path.read_text());cover=json.loads(cover_path.read_text());base=json.loads(base_path.read_text())
        node=cover['nodes'][663]
        if (node['id'],node['parent'],node['depth'],node['model_sha256'],node['local'],node['physical'])!=(
                663,[52,1],4,MODEL_SHA,LOCAL,PHYSICAL):raise ValueError('Source cover node changed.')
        if old['node_id']!=663 or old['local']!=LOCAL or old['physical']!=PHYSICAL or old['baseline_sha256']!=BASELINE_SHA:
            raise ValueError('Source model mismatch.')
        c1,h,R=map(Q,(base['c'],base['h'],base['R']))
        if (Q(base['alpha']),h,R,base['N'])!=(Q(47,40),Q(1,32),Q(33,20),80) or c1<=0:
            raise ValueError('Fixed baseline family mismatch.')
        u0,u1,v0,v1=map(Q,PHYSICAL);x0,x1,y0,y1=map(Q,LOCAL)
        uc,vc=(u0+u1)/2,(v0+v1)/2;xc,yc=(x0+x1)/2,(y0+y1)/2
        du,dv=(u1-u0)/2,(v1-v0)/2
        if (u0,u1,v0,v1)!=(h*(52+x0),h*(52+x1),h*(1+y0),h*(1+y1)) or not (u0>0 and v0>0 and u0+v0>R):
            raise ValueError('Physical geometry mismatch.')
        psi=old['Psi'];es=Q(old['AS']['remainder_upper']);muv=Q(psi['Muv_upper']);old_D=Q(psi['one_sided_lower_deduction'])
        if psi['one_sided_only'] is not True or psi['unscaled_by_c'] is not True or psi['Bernstein_degree']!=[3,3] or \
                es<0 or muv<0 or old_D!=muv*du*dv or Q(psi['du'])!=du or Q(psi['dv'])!=dv:
            raise ValueError('Source signed model metadata/arithmetic mismatch.')
        if len(psi['mixed_derivatives'])!=2:raise ValueError('Expected two stored Hab bounds.')
        aus=[]
        for br,(a,b) in zip(psi['mixed_derivatives'],((u0,v0),(v0,u0))):
            lo,hi,au=map(Q,(br['lower'],br['upper'],br['absolute_upper']))
            if (Q(br['a']),Q(br['b']),br['kind'])!=(a,b,'Hab') or lo>hi or au<max(abs(lo),abs(hi)):
                raise ValueError('Stored Hab absolute enclosure invalid.')
            aus.append(au)
        if sum(aus,Q(0))!=muv:raise ValueError('Stored Muv sum mismatch.')
        if any(len(old[key]['bernstein'])!=16 for key in ('Psi','AS')) or len(old['rows'])!=16:
            raise ValueError('Expected 16 source coefficients/rows.')
        t1=[];as_lower=[];old_c=[];old_s=[];old_slack=[]
        for i in range(4):
            tr=[];sr=[];cr=[];ssr=[];slar=[]
            for j in range(4):
                k=4*i+j;pb,sb,row=psi['bernstein'][k],old['AS']['bernstein'][k],old['rows'][k]
                if pb['index']!=[i,j] or sb['index']!=[i,j] or row['bernstein_index']!=[i,j] or \
                        Q(pb['lower'])>Q(pb['upper']) or Q(sb['lower'])>Q(sb['upper']):raise ValueError('Source indexing/interval invalid.')
                tr.append(Q(pb['lower']));sr.append(Q(sb['lower']));cr.append(Q(row['c_column']))
                ssr.append(Q(row['s_column']));slar.append(Q(row['baseline_slack']))
            t1.append(tr);as_lower.append(sr);old_c.append(cr);old_s.append(ssr);old_slack.append(slar)
        guard()
        rt,rs,roc,ros,rol=[tensor_quadrants(matrix) for matrix in (t1,as_lower,old_c,old_s,old_slack)]
        proof={'node_id':663,'source_model_sha256':MODEL_SHA,'source_cover_sha256':COVER_SHA,'baseline_sha256':BASELINE_SHA,
            'source_tree_depth':4,'source_tree_changed':False,'proof_partition':'four_center_sign_quadrants',
            'physical':PHYSICAL,'local':LOCAL,'physical_center':[str(uc),str(vc)],'local_center':[str(xc),str(yc)],
            'Muv_upper':str(muv),'old_uniform_D':str(old_D),'AS_remainder_upper':str(es),
            'source_T1_lower':strings(t1),'source_AS_lower':strings(as_lower),'pieces':[]}
        pu=((u0,uc),(uc,u1));pv=((v0,vc),(vc,v1));lx=((x0,xc),(xc,x1));ly=((y0,yc),(yc,y1))
        all_slacks=[];all_gains=[]
        for piece_id,(iu,iv) in enumerate(((0,0),(1,0),(0,1),(1,1))):
            guard();result['pieces_started']+=1;result['in_flight']={'piece_id':piece_id,'stage':'rows'};save()
            a,b=pu[iu];c,d=pv[iv];same_sign=iu==iv
            bilinear2=[[Q(0) if same_sign else muv*(u-uc)*(v-vc) for v in (c,d)] for u in (a,b)]
            correction=elevate_bilinear(bilinear2)
            piece={'piece_id':piece_id,'u_sign':'-' if iu==0 else '+','v_sign':'-' if iv==0 else '+',
                'proof_piece':True,'source_tree_depth':4,'physical':list(map(str,(a,b,c,d))),
                'local':list(map(str,(*lx[iu],*ly[iv]))),'same_sign':same_sign,
                'restricted_T1_lower':strings(rt[(iu,iv)]),'restricted_AS_lower':strings(rs[(iu,iv)]),
                'bilinear_Bernstein_degree1':strings(bilinear2),'bilinear_Bernstein_degree3':strings(correction),'rows':[]}
            for i in range(4):
                for j in range(4):
                    guard()
                    if result['rows_started']>=64:raise ValueError('New row budget exceeded.')
                    result['rows_started']+=1;t,s=rt[(iu,iv)][i][j],rs[(iu,iv)][i][j];q=correction[i][j]
                    cc,ss=t+q,s-es;slack=c1*cc+ss
                    oc,os,old_sl=roc[(iu,iv)][i][j],ros[(iu,iv)][i][j],rol[(iu,iv)][i][j]
                    if oc!=t-old_D or os!=ss or old_sl!=c1*oc+os:raise ValueError('Restricted old row relation mismatch.')
                    gain=slack-old_sl
                    if gain!=c1*(old_D+q) or gain<0:raise ValueError('Signed piece lower model failed dominance check.')
                    piece['rows'].append({'bernstein_index':[i,j],'c_column':str(cc),'s_column':str(ss),
                        'baseline_slack':str(slack),'restricted_old_c_column':str(oc),'restricted_old_s_column':str(os),
                        'restricted_old_slack':str(old_sl),'gain':str(gain)})
                    all_slacks.append((slack,piece_id,[i,j]));all_gains.append(gain);result['rows_completed']+=1
            piece['accepted']=all(Q(row['baseline_slack'])>0 for row in piece['rows'])
            proof['pieces'].append(piece);result['pieces_completed']+=1;result['in_flight']=None;save()
        guard();out=ROOT/'quadrants.json';dump(out,proof)
        worst=min(all_slacks,key=lambda item:item[0])
        result.update(status='completed',accepted=all(slack>0 for slack,_,_ in all_slacks),
            minimum_slack=str(worst[0]),minimum_piece=worst[1],minimum_index=worst[2],
            minimum_gain=str(min(all_gains)),quadrants_sha256=sha256(out.read_bytes()).hexdigest())
    except Exception as exc:
        result.update(status='failed',accepted=False,reason=f'{type(exc).__name__}: {exc}')
        if proof is not None:dump(ROOT/'quadrants.partial.json',proof)
    result['source_sha256']=sha256(Path(__file__).read_bytes()).hexdigest();save()
    print(json.dumps({'status':result['status'],'accepted':result.get('accepted',False),
        'pieces_completed':result['pieces_completed'],'rows_completed':result['rows_completed'],
        'elapsed_seconds':result['elapsed_seconds']}));return 0


if __name__=='__main__':sys.exit(main())
