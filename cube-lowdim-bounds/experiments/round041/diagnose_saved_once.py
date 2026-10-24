#!/usr/bin/env python3
"""One saved-model Fraction pass, at most 16 rows/60s. No mathematical library imports."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import signal
from time import monotonic

ROOT=Path(__file__).resolve().parent
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'
COVER_SHA='e78b2bf4e6f5e00974a180cb85dfa9044ed42939930c4add48078e6c4064c8e7'


def main():
    output=ROOT/'failure_breakdown.json'
    if output.exists():raise SystemExit('Existing diagnostic; no repeated pass.')
    start=monotonic();result={'status':'started','rows_started':0,'rows_completed':0,'maximum_rows':16,
        'seconds_hard':60,'new_AS':0,'new_Psi':0,'new_integrals':0,'new_points':0,'LP_calls':0,'rows':[]}
    def dump():
        result['elapsed_seconds']=monotonic()-start
        tmp=output.with_suffix('.json.tmp');tmp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');tmp.replace(output)
    def stop(signum,frame):raise TimeoutError('60s diagnostic deadline')
    def guard():
        if monotonic()-start>=60:raise TimeoutError('60s diagnostic deadline')
    def pair(value):return {'exact':str(value),'approx':float(value)}
    signal.signal(signal.SIGALRM,stop);signal.setitimer(signal.ITIMER_REAL,60)
    dump()
    try:
        cover_path=ROOT/'leaf_cover.json';model_path=ROOT/'models.local/node_00663.json'
        base_path=ROOT.parent/'round01/scaled_certificate.json'
        if sha256(cover_path.read_bytes()).hexdigest()!=COVER_SHA:raise ValueError('Frozen cover hash mismatch.')
        if sha256(base_path.read_bytes()).hexdigest()!=BASELINE_SHA:raise ValueError('Baseline hash mismatch.')
        cover=json.loads(cover_path.read_text());node=cover['nodes'][663]
        model_sha=sha256(model_path.read_bytes()).hexdigest()
        if node['id']!=663 or node['model_sha256']!=model_sha:raise ValueError('Model provenance mismatch.')
        model=json.loads(model_path.read_text());base=json.loads(base_path.read_text())
        if model['baseline_sha256']!=BASELINE_SHA or model['physical']!=node['physical'] or model['local']!=node['local']:
            raise ValueError('Node/model identity mismatch.')
        c1,h,R=map(Q,(base['c'],base['h'],base['R']))
        u0,u1,v0,v1=map(Q,model['physical']);x0,x1,y0,y1=map(Q,model['local'])
        if (u0,u1,v0,v1)!=(h*(52+x0),h*(52+x1),h*(1+y0),h*(1+y1)):
            raise ValueError('Local/physical mismatch.')
        if len(model['rows'])!=16 or len(model['AS']['bernstein'])!=16 or len(model['Psi']['bernstein'])!=16:
            raise ValueError('Expected exactly 16 saved rows.')
        es,ep=Q(model['AS']['remainder_upper']),Q(model['Psi']['one_sided_lower_deduction'])
        if es<0 or ep<0 or c1<=0:raise ValueError('Negative remainder or coefficient.')
        weighted_ep=c1*ep;total_e=weighted_ep+es
        result.update(inputs={'model_sha256':model_sha,'cover_sha256':COVER_SHA,'baseline_sha256':BASELINE_SHA},
            geometry={'node_id':663,'depth':model['depth'],'local':model['local'],'physical':model['physical'],
                'r_lower':str(u0+v0),'r_upper':str(u1+v1),'strict_exterior':u0+v0>R,
                'u_crosses_R':u0<R<u1,'v_crosses_R':v0<R<v1,
                'u_entirely_below_R':u1<R,'v_entirely_below_R':v1<R,'away_axes':u0>0 and v0>0},
            deductions_and_AS_remainder={'AS':pair(es),'Psi_one_sided':pair(ep),'weighted_Psi_one_sided':pair(weighted_ep),'total':pair(total_e)})
        psi=model['Psi'];du,dv=(u1-u0)/2,(v1-v0)/2
        if Q(psi['du'])!=du or Q(psi['dv'])!=dv or psi['one_sided_only'] is not True:
            raise ValueError('Saved halfwidth/single-sided label mismatch.')
        curvature=[];absolute_uppers=[];mixed_lowers=[];mixed_uppers=[]
        if len(psi['mixed_derivatives'])!=2:raise ValueError('Expected two saved Hab balls.')
        for dr,(a,b) in zip(psi['mixed_derivatives'],((u0,v0),(v0,u0))):
            lo,hi,au=map(Q,(dr['lower'],dr['upper'],dr['absolute_upper']))
            if (Q(dr['a']),Q(dr['b']),dr['kind'])!=(a,b,'Hab') or not lo<=hi or au<max(abs(lo),abs(hi)):
                raise ValueError('Saved mixed curvature interval invalid.')
            curvature.append({'a':str(a),'b':str(b),'lower':pair(lo),'upper':pair(hi),
                'ball_width':pair(hi-lo),'absolute_upper':pair(au)})
            absolute_uppers.append(au);mixed_lowers.append(lo);mixed_uppers.append(hi)
        muv=sum(absolute_uppers,Q(0))
        if Q(psi['Muv_upper'])!=muv or ep!=muv*du*dv:raise ValueError('Saved D arithmetic mismatch.')
        center={}
        for key in ('F0','Fu','Fv'):
            data=psi['center_values'][key];lo,hi=Q(data['lower']),Q(data['upper'])
            if lo>hi:raise ValueError('Reversed saved center interval.')
            center[key]={'lower':pair(lo),'upper':pair(hi),'ball_width':pair(hi-lo)}
        result.update(saved_center_T1_coefficients=center,saved_mixed_derivatives=curvature,
            mixed_curvature={'Muv_upper':pair(muv),'sum_Hab_lower':pair(sum(mixed_lowers,Q(0))),
                'sum_Hab_upper':pair(sum(mixed_uppers,Q(0))),
                'sum_ball_width':pair(sum(mixed_uppers,Q(0))-sum(mixed_lowers,Q(0))),
                'D_due_to_ball_width':pair((sum(mixed_uppers,Q(0))-sum(mixed_lowers,Q(0)))*du*dv),
                'du':str(du),'dv':str(dv)})
        aggregate=[]
        for ell in range(16):
            guard();result['rows_started']+=1
            pb,sb,row=model['Psi']['bernstein'][ell],model['AS']['bernstein'][ell],model['rows'][ell]
            index=[ell//4,ell%4]
            if pb['index']!=index or sb['index']!=index or row['bernstein_index']!=index:raise ValueError('Index mismatch.')
            pl,pu,sl,su=map(Q,(pb['lower'],pb['upper'],sb['lower'],sb['upper']))
            if pl>pu or sl>su:raise ValueError('Reversed saved interval.')
            lower,upper=c1*pl+sl,c1*pu+su
            after=lower-total_e
            if Q(row['c_column'])!=pl-ep or Q(row['s_column'])!=sl-es or Q(row['baseline_slack'])!=after:
                raise ValueError('Saved row arithmetic mismatch.')
            result['rows'].append({'index':index,'Psi_lower':str(pl),'Psi_upper':str(pu),
                'AS_lower':str(sl),'AS_upper':str(su),'Psi_ball_width':pair(pu-pl),'AS_ball_width':pair(su-sl),
                'combined_ball_width':pair(upper-lower),'before_remainder_lower':pair(lower),
                'before_remainder_upper':pair(upper),'after_remainder_lower':pair(after),
                'before_nonpositive':lower<=0,'after_nonpositive':after<=0})
            aggregate.append((index,lower,upper,after,upper-lower));result['rows_completed']+=1
        guard()
        worst=min(aggregate,key=lambda x:x[3]);best=max(aggregate,key=lambda x:x[1])
        result.update(status='completed',summary={'minimum_index':worst[0],
            'minimum_before_remainder':pair(worst[1]),'minimum_after_remainder':pair(worst[3]),
            'maximum_coefficient_lower':pair(best[1]),
            'coefficient_lower_spread':pair(best[1]-worst[1]),
            'coefficient_interval_hull_width':pair(max(x[2] for x in aggregate)-min(x[1] for x in aggregate)),
            'maximum_combined_ball_width':pair(max(x[4] for x in aggregate)),
            'before_nonpositive_rows':sum(x[1]<=0 for x in aggregate),
            'after_nonpositive_rows':sum(x[3]<=0 for x in aggregate)})
    except Exception as exc:result.update(status='failed',reason=f'{type(exc).__name__}: {exc}')
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        result['source_sha256']=sha256(Path(__file__).read_bytes()).hexdigest();dump()
    print(json.dumps({key:result[key] for key in ('status','rows_started','rows_completed','elapsed_seconds')}))


if __name__=='__main__':main()
