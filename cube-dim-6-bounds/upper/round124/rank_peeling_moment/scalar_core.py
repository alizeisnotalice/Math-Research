"""Exact scalar model and interval certificates; no matrix/sample computation."""
from fractions import Fraction as Q
from math import comb,factorial

def req(x,label):
    if not x:raise ArithmeticError(label)
def down(q,bits):
    q=Q(q);return Q((q.numerator*(1<<bits))//q.denominator,1<<bits)
def up(q,bits):
    q=Q(q);return Q(-((-q.numerator*(1<<bits))//q.denominator),1<<bits)
def phi(n,k):
    if k==0:return 1
    if k<=3:return comb(n,k) if n>=k else 0
    last=1<<(k-2)
    if n<=last:return 0
    product=n
    for j in range(k-1):product*=n-(1<<j)
    return (product+factorial(k)-1)//factorial(k)
def pmfs(p):
    return [Q(comb(32,n))*p**n*(1-p)**(32-n) for n in range(33)]
def displays(iv):return [format(float(Q(x)),'.17g') for x in iv]

class Core:
    def __init__(self,inp):
        from flint import arb,fmpq,ctx
        ctx.prec=512;self.arb=arb;self.fmpq=fmpq
        self.threshold=Q(inp['threshold_M_power6_exact'])
        req(self.threshold==Q(8,9)**6*24576**5/65000**6,'complete123 exact threshold identity')
        R0=Q(inp['full_fan_R0']);R1=Q(inp['R1'])
        req(R0==Q(inp['independent_fan_R0']) and R0<Q(27,2**21),'old exact full fan improves108 B0')
        req([c['unordered_count'] for c in inp['R3full']]==[480,1920,480,640,1440], 'full4960, never no-ones count')
        req(sum(c['unordered_count'] for c in inp['R3full'])==4960,'all full unordered triples')
        self.rational_coeff={0:R0,1:64*R1,
            2:4*sum(Q(c['unordered_count'])*Q(c['R']) for c in inp['R2']),
            3:8*sum(Q(c['unordered_count'])*Q(c['R']) for c in inp['R3full']),
            4:Q(comb(32,4)*5,192),6:Q(2*comb(32,6))}
        self.table=[[phi(n,k) for k in range(7)] for n in range(33)]
        for k in range(7):
            req(all(self.table[n][k]<=self.table[n+1][k] for n in range(32)),'Phi monotonicity')
        for n,row in enumerate(self.table):
            req(all(0<=v<=(comb(n,k) if n>=k else 0) for k,v in enumerate(row)),'legal integer count upper sanity')

    def A(self,q):
        q=Q(q);return self.arb(self.fmpq(q.numerator,q.denominator))
    def ends(self,x):
        vals=[]
        for z in (x.lower(),x.upper()):
            m,e=z.man_exp();vals.append(Q(int(m))*Q(2)**int(e))
        return vals
    def grid(self,x,bits):
        l,u=self.ends(x);return down(l,bits),up(u,bits)
    def coefficients(self):
        return {'pi_cubed_B_rational':{str(k):str(v) for k,v in self.rational_coeff.items()},
                'pi_cubed_B5':'binom(32,5)*sqrt(8/3)/8',
                'B0_is_exact_full_fan':True,'B4_to_B6_are_coarse_independent_subset_caps':True,
                'Phi_table':self.table,'complete_threshold_M_power6':str(self.threshold)}
    def caps(self,h):
        a=self.A(h);pi=self.arb.pi();z=1-(-a).exp();p=a.sqrt().erf()
        pL,pU=self.grid(p,192);req(0<pL<pU<1,'actual binomial probability')
        beta=[];lower=[]
        for k in range(7):
            coef=self.A(self.rational_coeff[k]) if k!=5 else self.A(Q(comb(32,5),8))*self.A(Q(8,3)).sqrt()
            b=coef/pi**3*(z**k if k else 1)
            lo,hi=self.grid(b,240);req(0<lo<=hi,'positive cap enclosure')
            lower.append(lo);beta.append(hi)
        return {'h':str(h),'p_interval':[str(pL),str(pU)],
                'b_lower':[str(x) for x in lower],'beta_upper':[str(x) for x in beta],
                'p_display':displays([pL,pU]),'cap_rounding_absolute_bits':240,
                'p_rounding_absolute_bits':192}
    def certify(self,caps,integers):
        den=1<<44;req(len(integers)==7 and sum(integers)==den and integers[0]>0 and min(integers)>=0,'exact rational simplex')
        x=[Q(v,den) for v in integers];beta=list(map(Q,caps['beta_upper']));bl=list(map(Q,caps['b_lower']))
        pL,pU=map(Q,caps['p_interval']);probsL,probsU=pmfs(pL),pmfs(pU)
        req(sum(probsL)==sum(probsU)==1,'both endpoint binomial laws normalized exactly')
        probmax=[]
        for n in range(33):
            values=[probsL[n],probsU[n]]
            critical=Q(n,32)
            if pL<=critical<=pU:values.append(Q(comb(32,n))*critical**n*(1-critical)**(32-n))
            probmax.append(max(values))
        features=[[Q(v)/b for v,b in zip(row,beta)] for row in self.table]
        roots=[]
        for n,row in enumerate(features):
            w=sum(v*y for v,y in zip(row,x));req(w>0,'W positive on every count state')
            vl,vu=self.grid(1/self.A(w).root(5),192)
            req(0<vl<=vu and vl**5<=1/w<=vu**5,'exact reciprocal fifth-root enclosure')
            roots.append({'n':n,'W':str(w),'v_interval':[str(vl),str(vu)]})
        Hlo=sum(probsU[n]*Q(roots[n]['v_interval'][0]) for n in range(33))
        Hup=sum(probsL[n]*Q(roots[n]['v_interval'][1]) for n in range(33))
        req(0<Hlo<=Hup,'positive33-term full binomial expectation')
        gradients=[]
        for k in range(7):
            gradients.append(-Q(1,5)*sum(probmax[n]*features[n][k]*Q(roots[n]['v_interval'][1])**6 for n in range(33)))
        dual=Q(6,5)*Hlo+min(gradients);gamma=min(a/b for a,b in zip(bl,beta))
        upper6=Hup**5;familylower6=gamma*max(Q(0),dual)**5
        ml,mu=self.grid(self.A(upper6).root(6),192)
        req(0<ml<=mu and ml**6<=upper6<=mu**6,'exact sixth root of certificate expression')
        status=('STRICT_MOMENT_THRESHOLD_PASS' if upper6<self.threshold else
                'FIXED_SINGLE_H_TABLE_RIGOROUSLY_OBSTRUCTED' if familylower6>self.threshold else
                'UNSEPARATED_FIXED_CANDIDATE')
        return {'status':status,'x_grid_denominator':str(den),'x_grid_numerators':integers,
                'a_rational':[str(v/b) for v,b in zip(x,beta)],'EdW_upper_exact':'1',
                'root_certificates':roots,'H_interval':[str(Hlo),str(Hup)],
                'gradient_lower':[str(v) for v in gradients],'support_plane_lower':str(dual),
                'gamma_exact_cap_rescaling_lower':str(gamma),
                'M_power6_upper':str(upper6),'threshold_minus_upper_exact':str(self.threshold-upper6),
                'fixed_cap_family_M_power6_lower':str(familylower6),
                'fixed_cap_family_lower_minus_threshold':str(familylower6-self.threshold),
                'M_upper_expression_interval':[str(ml),str(mu)],'M_upper_expression_display':displays([ml,mu]),
                'M_power6_upper_over_threshold_display':format(float(upper6/self.threshold),'.17g'),
                'full_I_upper_power5_from123':str(Q(3,16)**6*upper6/384**5),
                'actual_M_lower_bound':False,'global_g0_proved':False,'new_C6_bound':False}
