"""Strict true-corner norm brackets using r028 exact reduction and directed Gamma integrals."""
from pathlib import Path
import multiprocessing as mp,json,sys,time,math
from functools import lru_cache
from certify82 import D,F,C,Z,O,T,exact,add,mul,div,logpoint,exp_iv,gamma

def expsigned(x):
    if x[0]>=0:return exp_iv(x)
    if x[1]<=0:return div(exact(1),exp_iv((x[1].copy_negate(),x[0].copy_negate())))
    return F.divide(O,exp_iv(exact(x[0].copy_negate()))[1]),exp_iv(exact(x[1]))[1]

def run(z):
    n=z['N'];A=z['A'];b=D(z['candidate']);side=z['side'];started=time.monotonic()
    if n<2 or n&(n-1) or A not in [2,3] or b<=0 or side not in ['lower','upper']:raise ValueError('unsupported parameters')
    lnN=mul(exact(n.bit_length()-1),logpoint(T))
    @lru_cache(None)
    def node(t):
        l=logpoint(t);u=logpoint(F.add(O,l[0]))[0],logpoint(C.add(O,l[1]))[1]
        inner=add(exact(1),div(l,exact(n)));vsmall=logpoint(inner[0])[0],logpoint(inner[1])[1]
        v=add(mul(exact(n-2),u),add(lnN,vsmall));s=mul(exact(b),exp_iv(div(v,exact(n-1))))
        qq,ff=gamma(n,s)
        sp=div(mul(s,add(exact(n-1),l)),mul(exact(t),mul(add(exact(1),l),add(exact(n),l))))
        exponent=(F.subtract(v[0],C.subtract(s[1],b)),C.subtract(v[1],F.subtract(s[0],b)))
        return dict(s=s,Q=qq,pdf=ff,sprime=sp,cost=expsigned(exponent))
    baseQ=gamma(n,exact(b))[0];modepdf=gamma(n,exact(n-1))[1]
    kmode=None
    if b<D(n-1):
        ratio=div(exact(n-1),exact(b));lg=logpoint(ratio[0])[0],logpoint(ratio[1])[1]
        v=add(mul(exact(n-1),lg),(F.subtract(b,D(n-1)),C.subtract(b,D(n-1))));kmode=expsigned(v)[1]
    for power in range(5,14):
        M=2**power;h=F.divide(D(A-1),D(M));il=ih=cl=ch=Z
        for i in range(M):
            a=F.add(O,F.multiply(D(i),h));bb=F.add(a,h);mid=F.divide(F.add(a,bb),T);na,nb,nm=node(a),node(bb),node(mid)
            fl=min(na['pdf'][0],nb['pdf'][0]);fh=max(na['pdf'][1],nb['pdf'][1])
            crosses=na['s'][0]<=n-1<=nb['s'][1]
            if crosses:fh=max(fh,modepdf[1])
            dl=F.multiply(fl,nb['sprime'][0]);dh=C.multiply(fh,na['sprime'][1]);err=C.divide(C.multiply(C.multiply(h,h),C.subtract(dh,dl)),D(8))
            il=F.add(il,F.subtract(F.multiply(h,nm['Q'][0]),err));ih=C.add(ih,C.add(C.multiply(h,nm['Q'][1]),err))
            kl=min(na['cost'][0],nb['cost'][0]);kh=max(na['cost'][1],nb['cost'][1])
            if crosses and kmode is not None:kh=max(kh,kmode)
            cl=F.add(cl,F.multiply(h,kl));ch=C.add(ch,C.multiply(h,kh))
        ml=F.add(baseQ[0],il);mh=C.add(baseQ[1],ih)
        if (side=='upper' and mh<O) or (side=='lower' and ml>O):break
    else:raise ArithmeticError(('mean sign not certified',z,str(ml),str(mh)))
    zl=F.add(O,cl);zh=C.add(O,ch)
    loglo=F.subtract(b,logpoint(zh)[1]);loghi=C.subtract(b,logpoint(zl)[0])
    return dict(z,partition_cells=M,evaluated_nodes=node.cache_info().currsize,mean_interval=[str(ml),str(mh)],scaled_modular_interval=[str(zl),str(zh)],log_norm_at_candidate_interval=[str(loglo),str(loghi)],
      proof='directed positive-series Gamma tails; midpoint derivative-range bound for mean; unimodal positive cost rectangles',pid=__import__('os').getpid(),seconds=time.monotonic()-started)

def main():
    src=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True);rows=json.loads(src.read_text())['rows'];jobs=[]
    for r in rows:
        if r['N']>256:continue
        for side,v in [('lower',(math.floor(100*r['activation'])-1)/100),('upper',(math.ceil(100*r['activation'])+1)/100)]:jobs.append(dict(N=r['N'],A=r['A'],candidate=str(D(str(v))),side=side))
    started=time.monotonic()
    with mp.get_context('spawn').Pool(12) as pool:cert=list(pool.map(run,jobs))
    brackets=[]
    for n in [2,4,8,16,32,64,128,256]:
        for A in [2,3]:
            lo=next(r for r in cert if r['N']==n and r['A']==A and r['side']=='lower');hi=next(r for r in cert if r['N']==n and r['A']==A and r['side']=='upper')
            brackets.append(dict(N=n,A=A,activation_interval=[lo['candidate'],hi['candidate']],log_true_norm_interval=[lo['log_norm_at_candidate_interval'][0],hi['log_norm_at_candidate_interval'][1]]))
    result=dict(pass_all=True,workers=12,actual_processes=len({r['pid'] for r in cert}),precision=70,seconds=time.monotonic()-started,rows=cert,brackets=brackets)
    (out/'certificates83.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(pass_all=True,brackets=len(brackets),workers=result['actual_processes'],cells=sum(r['partition_cells'] for r in cert),max_cells=max(r['partition_cells'] for r in cert),seconds=result['seconds'],N64cap3=next(r for r in brackets if r['N']==64 and r['A']==3))))
if __name__=='__main__':main()
