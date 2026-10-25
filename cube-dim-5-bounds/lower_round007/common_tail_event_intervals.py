#!/usr/bin/env python3
"""Local exact event-interval certificate for the common-last-four parameter.

No displayed decimal is used in the construction.  The finite endpoint event
list contains every affine endpoint alternative before box dominance removal.
Within an event interval all endpoint orders and cell masks are fixed; hence
the exact union ratio is a polynomial of degree at most four in delta.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, prod
from pathlib import Path
import hashlib, json, sys, time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path[:0]=[str(ROOT/'lower'),str(ROOT/'lower_round004')]
from certify_rectangular import field
from common_tail_batch import root_value, evaluate

WINDOW=Q(1,16); MAX_SIDE=4; MAX_SUBDIVISIONS=32
PER_EVAL_SECONDS=120

def lift(F,q): return F(q.numerator)/q.denominator
def enc(x): return {f'{i//5},{i%5}':str(c) for i,c in enumerate(x.a) if c}
def dec(F,d):
    a=[Q(0)]*25
    for key,val in d.items():
        i,j=map(int,key.split(','));a[5*i+j]=Q(val)
    return F(a=a)

def endpoint_forms(F,u,v):
    """All endpoint alternatives b+m*delta, pre-min/max and pre-dominance."""
    a1=1+u**3/2; base=(1+u)/2
    forms={(F(0),Q(0)),(a1/2,Q(0))}
    for ls in product(range(1,5),repeat=5):
        r=root_value(prod(ls),u,v)/2
        # first axis is fixed; include its selectable raw endpoints.
        l=ls[0]
        if l%2: forms.add((r-(l-1)*a1/2,Q(0)))
        else: forms.add((l*a1/2-r,Q(0)))
        # Every tail axis has the same possible formulas, but r depends on ls.
        for l in ls[1:]:
            forms.add((base/2,Q(1,2)))
            if l%2: forms.add((r-(l-1)*base/2,Q(-(l-1),2)))
            else: forms.add((l*base/2-r,Q(l,2)))
    return tuple(forms)

def event_values(F,u,v,bounds):
    forms=endpoint_forms(F,u,v)
    low,high=lift(F,-WINDOW),lift(F,WINDOW)
    ev={low.a:low, high.a:high, F(0).a:F(0)}
    for (b,m),(c,n) in combinations(forms,2):
        if m==n: continue
        x=(c-b)/(m-n)  # slope denominator is rational
        if low < x and x < high: ev[x.a]=x
    return forms,tuple(sorted(ev.values()))

def rational_nodes(left,right,bounds):
    # Strict rational enclosure from isolating intervals, then six distinct
    # interior rational points.  Their ordering is checked in Q(u,v).
    _,lh=bounds(left.a); rl,_=bounds(right.a)
    if not lh < rl: raise ArithmeticError('event isolating intervals do not separate')
    nodes=[lh+(rl-lh)*Q(k,7) for k in range(1,7)]
    return nodes

def interpolate(F,xs,ys):
    out=[F(0)]*5
    for i,(x,y) in enumerate(zip(xs[:5],ys[:5])):
        poly=[F(1)]; den=Q(1)
        for j,z in enumerate(xs[:5]):
            if i==j:continue
            nxt=[F(0)]*(len(poly)+1)
            for k,c in enumerate(poly):
                nxt[k] += c*(-z); nxt[k+1] += c
            poly=nxt; den*=x-z
        for k,c in enumerate(poly): out[k]+=y*c/den
    return out

def poly_value(coeff,x):
    out=coeff[-1]
    for c in reversed(coeff[:-1]): out=out*x+c
    return out

def bernstein(coeff,degree,left,right):
    """Bernstein coefficients of sum coeff[i] x^i on [left,right]."""
    h=right-left; q=[coeff[0]*0 for _ in range(degree+1)]
    for i,c in enumerate(coeff):
        for k in range(i+1): q[k]+=c*Q(comb(i,k))*left**(i-k)*h**k
    return [sum((q[i]*Q(comb(k,i),comb(degree,i)) for i in range(k+1)),q[0]*0) for k in range(degree+1)]

def nonpositive_on(coeff,left,right,F):
    # Bounded breadth-first subdivision; Bernstein convex-hull property.
    todo=[(left,right)]; used=0
    while todo:
        lo,hi=todo.pop(); bs=bernstein(coeff,4,lo,hi)
        if all(b <= F(0) for b in bs): continue
        if used>=MAX_SUBDIVISIONS: return False,used,len(todo)+1
        mid=(lo+hi)/2; todo.extend(((lo,mid),(mid,hi))); used+=1
    return True,used,0

def choose_intervals(events,F):
    zero=F(0)
    left=[x for x in events if x<zero]
    right=[x for x in events if zero<x]
    # Splitting at C is allowed even when it is not an event.  It gives exactly
    # the requested four nearest one-sided intervals, not a global window proof.
    # Keep only the requested nearest four *intervals*.  Adding the window
    # endpoint when four events already exist would silently span older events.
    lb=(left[-MAX_SIDE:] if len(left)>=MAX_SIDE else [lift(F,-WINDOW)]+left)+[zero]
    rb=[zero]+(right[:MAX_SIDE] if len(right)>=MAX_SIDE else right+[lift(F,WINDOW)])
    def distinct_pairs(x):
        x=sorted({z.a:z for z in x}.values())
        return list(zip(x,x[1:]))
    return [('left',a,b) for a,b in distinct_pairs(lb)]+[('right',a,b) for a,b in distinct_pairs(rb)]

def analyze_interval(F,u,v,bounds,c_upper,side,left,right):
    xs=rational_nodes(left,right,bounds)
    if not all(left < lift(F,x) and lift(F,x) < right for x in xs): raise AssertionError('node outside event interval')
    vals=[]; evidence=[]
    for x in xs:
        row=evaluate(F,u,v,bounds,x,time.monotonic()+PER_EVAL_SECONDS)
        vals.append(dec(F,row['ratio_coefficients'])); evidence.append(row)
    coeff=interpolate(F,xs,vals)
    assert poly_value(coeff,xs[5])==vals[5]
    comparison=list(coeff); comparison[0]-=c_upper
    excluded,subs,unresolved=nonpositive_on(comparison,left,right,F)
    # Derivative sign is auxiliary: it explains local monotonicity only.
    derivative=[coeff[i]*i for i in range(1,5)]
    dpoly=derivative+[F(0)]
    neg,_,_=nonpositive_on(dpoly,left,right,F)
    pos,_,_=nonpositive_on([-z for z in dpoly],left,right,F)
    return {'side':side,'left_delta_coefficients':enc(left),'right_delta_coefficients':enc(right),
      'sample_deltas':[str(x) for x in xs], 'sample_union_certificates':evidence,
      'ratio_polynomial_coefficients':[enc(c) for c in coeff],
      'degree_at_most_four_validated_at_six_nodes':True,
      'no_strict_improvement_vs_exact_C_by_nonpositive_bernstein':excluded,'bernstein_subdivisions':subs,'unresolved_cells':unresolved,
      'derivative_nonpositive_by_bernstein':neg,'derivative_nonnegative_by_bernstein':pos}

def interval_metadata(F,u,v,bounds,index):
    forms,events=event_values(F,u,v,bounds)
    intervals=choose_intervals(events,F)
    side,left,right=intervals[index]
    return forms,events,side,left,right,rational_nodes(left,right,bounds)

def main():
    started=time.monotonic(); F,u,v,bounds,_,_=field()
    forms,events=event_values(F,u,v,bounds)
    intervals=choose_intervals(events,F)
    lock=json.loads((ROOT/'lower_round003/lock_c_a2_equals_a3.json').read_text())
    c_source=json.loads((ROOT/'lower_round003/a2_equals_a3_results.json').read_text())
    c_exact=dec(F,c_source['coefficients'])
    # One interval per invocation keeps every finite calculation visible and
    # recoverable under its own 120-second checkpoints.
    if len(sys.argv)==3 and sys.argv[1]=='--index':
        index=int(sys.argv[2]); side,left,right=intervals[index]
        row=analyze_interval(F,u,v,bounds,c_exact,side,left,right)
        out={'status':'completed_single_event_interval','index':index,'raw_affine_endpoint_form_count':len(forms),
             'event_count_including_window_and_C':len(events),'row':row,'elapsed_seconds':time.monotonic()-started,
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (HERE/f'event_interval_{index}_results.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(out,indent=2)); return
    # Sampling is separately resumable because six dual-union calls can exceed
    # an interactive terminal slice despite each individual call staying under
    # its registered 120-second limit.
    if len(sys.argv)==4 and sys.argv[1]=='--sample':
        index,k=int(sys.argv[2]),int(sys.argv[3]); forms,events,side,left,right,xs=interval_metadata(F,u,v,bounds,index)
        row=evaluate(F,u,v,bounds,xs[k],time.monotonic()+PER_EVAL_SECONDS)
        out={'status':'completed_exact_sample','index':index,'sample_index':k,'side':side,
             'left_delta_coefficients':enc(left),'right_delta_coefficients':enc(right),'delta':str(xs[k]),
             'union_certificate':row,'elapsed_seconds':time.monotonic()-started,
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (HERE/f'corrected_event_interval_{index}_sample_{k}.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(out,indent=2)); return
    if len(sys.argv)==3 and sys.argv[1]=='--finalize':
        index=int(sys.argv[2]); forms,events,side,left,right,xs=interval_metadata(F,u,v,bounds,index)
        samples=[json.loads((HERE/f'corrected_event_interval_{index}_sample_{k}.json').read_text()) for k in range(6)]
        if any(Q(x['delta'])!=xs[k] for k,x in enumerate(samples)): raise AssertionError('sample provenance mismatch')
        vals=[dec(F,x['union_certificate']['ratio_coefficients']) for x in samples]
        coeff=interpolate(F,xs,vals); assert poly_value(coeff,xs[5])==vals[5]
        comparison=list(coeff); comparison[0]-=c_exact
        excluded,subs,unresolved=nonpositive_on(comparison,left,right,F)
        dpoly=[coeff[i]*i for i in range(1,5)]+[F(0)]
        neg,_,_=nonpositive_on(dpoly,left,right,F); pos,_,_=nonpositive_on([-z for z in dpoly],left,right,F)
        row={'side':side,'left_delta_coefficients':enc(left),'right_delta_coefficients':enc(right),'sample_deltas':[str(x) for x in xs],
             'sample_union_certificates':[x['union_certificate'] for x in samples], 'ratio_polynomial_coefficients':[enc(c) for c in coeff],
             'degree_at_most_four_validated_at_six_nodes':True,'no_strict_improvement_vs_exact_C_by_nonpositive_bernstein':excluded,
             'bernstein_subdivisions':subs,'unresolved_cells':unresolved,'derivative_nonpositive_by_bernstein':neg,'derivative_nonnegative_by_bernstein':pos}
        out={'status':'finalized_event_interval_exact_C','index':index,'raw_affine_endpoint_form_count':len(forms),'event_count_including_window_and_C':len(events),'row':row,
             'elapsed_seconds':time.monotonic()-started,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        (HERE/f'exactC_event_interval_{index}_results.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(out,indent=2)); return
    if len(sys.argv)!=1: raise SystemExit('usage: common_tail_event_intervals.py [--index 0..7 | --sample i k | --finalize i]')
    rows=[]
    for side,left,right in intervals: rows.append(analyze_interval(F,u,v,bounds,c_upper,side,left,right))
    out={'status':'completed_local_event_interval_audit','parameter':'delta where tail=(1+u)/2+delta',
      'event_window':str(WINDOW),'raw_affine_endpoint_form_count':len(forms),'event_count_including_window_and_C':len(events),
      'selected_interval_count':len(rows),'selection':'up to four nearest intervals on each side of C=delta0',
      'C_exact_coefficients_source':'lower_round003/a2_equals_a3_results.json', 'intervals':rows,
      'all_selected_intervals_excluded':all(r['no_strict_improvement_vs_exact_C_by_nonpositive_bernstein'] for r in rows),
      'elapsed_seconds':time.monotonic()-started,
      'scope':'Only the listed local event intervals. This is not a global common-tail optimization or a C5 transfer proof.',
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'C_lock_sha256':hashlib.sha256((ROOT/'lower_round003/lock_c_a2_equals_a3.json').read_bytes()).hexdigest()}
    (HERE/'common_tail_event_intervals_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
