"""Exact positive multinomial sums and saved reciprocal-root certificates."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from functools import lru_cache
from math import comb,factorial,prod
import gzip,hashlib,json,os

def req(x,s):
    if not x:raise ArithmeticError(s)
def floor_scaled(q,S):q=Q(q);return (q.numerator*S)//q.denominator
def ceil_scaled(q,S):q=Q(q);return -((-q.numerator*S)//q.denominator)

class Joint:
    def __init__(self,inp,candidate):
        from flint import arb,ctx
        ctx.prec=192;self.arb=arb;self.S=1<<192;self.invS=arb(1)/self.S;self.D=self.S**33;self.S6=self.S**6
        self.phi=inp['Phi_table'];self.beta0=Q(inp['beta0']);self.beta=[[Q(z) for z in row]for row in inp['beta_by_h']]
        grid=int(candidate['denominator']);nums=candidate['numerators']
        req(grid==1<<44 and len(nums)==49 and min(nums)>=0 and nums[0]>0 and sum(nums)==grid,'positive exact49 simplex')
        self.x=[Q(v,grid) for v in nums];self.a0=self.x[0]/self.beta0
        self.a=[[self.x[1+6*j+k]/self.beta[j][k] for k in range(6)]for j in range(8)]
        self.common=[floor_scaled(self.a0,self.S),ceil_scaled(self.a0,self.S)]
        req(self.common[0]>0,'positive count0 lower bound')
        self.table=[]
        for j in range(8):
            tab=[]
            for n in range(33):
                value=sum(self.a[j][k]*self.phi[n][k+1] for k in range(6))
                tab.append([floor_scaled(value,self.S),ceil_scaled(value,self.S)])
            self.table.append(tab)
        self.fact=[factorial(k)for k in range(33)]
        self.p=[];self.q=[];self.pows=[]
        for side in range(2):
            ps=[]
            for pair in inp['p_intervals']:
                v=Q(pair[side])*self.S;req(v.denominator==1,'saved192-grid cumulative p');ps.append(int(v))
            qs=[ps[0]]+[ps[j]-ps[j-1]for j in range(1,8)]+[self.S-ps[-1]]
            req(min(qs)>0 and sum(qs)==self.S,'one honest whole endpoint law')
            self.p.append(ps);self.q.append(qs);self.pows.append([[z**n for n in range(33)]for z in qs])
        req(all(a<=b for a,b in zip(self.p[0],self.p[1])),'componentwise CDF coupling')
        self.root_requests=0
        @lru_cache(maxsize=65536)
        def root(wl,wu):
            req(0<wl<=wu,'positive integer W interval')
            hull=(arb(wl)*self.invS).union(arb(wu)*self.invS)
            v=1/hull.root(5)
            out=[]
            for side,z in enumerate((v.lower(),v.upper())):
                m,e=z.man_exp();m,e=int(m),int(e)+192
                if e>=0:value=m<<e
                elif side==0:value=m>>(-e)
                else:value=(m+(1<<(-e))-1)>>(-e)
                out.append(value)
            l,u=out;req(0<l<=u and l**5*wu<=self.S6<=u**5*wl,'exact state reciprocal fifth-root certificate')
            return l,u
        self.root=root
    def tables(self):
        return {'scale_power2':192,'common_W_units':list(map(str,self.common)),
                'group_W_units':[[list(map(str,z))for z in row]for row in self.table],
                'p_units':[[str(z)for z in row]for row in self.p],
                'category_units':[[str(z)for z in row]for row in self.q],
                'probability_denominator_power2':6144,'weighted_sum_denominator_power2':6336,
                'a0':str(self.a0),'a_by_h':[[str(z)for z in row]for row in self.a]}
    def w(self,state):
        return tuple(self.common[side]+sum(self.table[j][n][side]for j,n in enumerate(state))for side in (0,1))
    def probnum(self,state,side):
        cats=[state[0]]+[state[j]-state[j-1]for j in range(1,8)]+[32-state[-1]]
        coef=self.fact[32]//prod(self.fact[c]for c in cats)
        return coef*prod(self.pows[side][j][c]for j,c in enumerate(cats))
    def certify(self,rootdir,expected_hash,atomic,progress):
        rootdir=Path(rootdir);chunks=rootdir/'cert_chunks';chunks.mkdir(exist_ok=False)
        atomic('integer_tables.json',self.tables())
        headers=[];globalhash=hashlib.sha256();chunkhash=hashlib.sha256();records=[]
        totalL=totalU=0;partialL=partialU=0;countK=[0]*17;partK=[0]*17
        totalMass=[0,0];partMass=[0,0]
        first=None;chunkfirst=None;last=None;index=0;start=0
        def flush():
            nonlocal records,partialL,partialU,partK,partMass,chunkhash,start,chunkfirst
            if not records:return
            k=len(headers);name='roots_'+str(k).zfill(3)+'.json.gz';path=chunks/name;tmp=path.with_name(path.name+'.own-tmp')
            payload=json.dumps(records,separators=(',',':')).encode('ascii')
            with tmp.open('wb')as raw:
                with gzip.GzipFile(filename='',mode='wb',fileobj=raw,compresslevel=6,mtime=0)as f:f.write(payload)
                raw.flush();os.fsync(raw.fileno())
            os.replace(tmp,path);req(path.stat().st_size<100_000_000,'compressed file project size')
            header={'chunk':k,'root_file':'cert_chunks/'+name,'root_file_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                    'start':start,'end_exclusive':index,'states':len(records),'first':list(chunkfirst),'last':list(last),
                    'state_bytes_sha256':chunkhash.hexdigest(),'prefix_count_by_K':partK,
                    'lower_at_pU_numerator':str(partialL),'upper_at_pL_numerator':str(partialU),
                    'prefix_probability_numerator_by_law':list(map(str,partMass))}
            atomic('cert_chunks/header_'+str(k).zfill(3)+'.json',header);headers.append(header)
            atomic('certificate_checkpoint.json',{'phase':'single_candidate_integer_prefix','completed_states':index,
                'chunks':headers,'lower_numerator':str(totalL),'upper_numerator':str(totalU),'prefix_count_by_K':countK})
            progress(index,len(headers))
            start=index;records=[];partialL=partialU=0;partK=[0]*17;partMass=[0,0];chunkhash=hashlib.sha256();chunkfirst=None
        for state in combinations_with_replacement(range(17),8):
            if first is None:first=state
            if chunkfirst is None:chunkfirst=state
            last=state;b=bytes(state);globalhash.update(b);chunkhash.update(b)
            wl,wu=self.w(state);vl,vu=self.root(wl,wu);self.root_requests+=1
            pl,pu=self.probnum(state,0),self.probnum(state,1)
            low=pu*vl;upp=pl*vu
            totalMass[0]+=pl;totalMass[1]+=pu;partMass[0]+=pl;partMass[1]+=pu
            totalL+=low;totalU+=upp;partialL+=low;partialU+=upp;countK[state[-1]]+=1;partK[state[-1]]+=1
            records.append([str(vl),str(vu)]);index+=1
            if len(records)==8192:flush()
        flush()
        req(index==735471 and countK==[comb(k+7,7)for k in range(17)],'complete735471 prefix')
        req(globalhash.hexdigest()==expected_hash,'float and exact complete state order match')
        tail=[];tailU=0;tailMass=[0,0];p=self.p[0][-1]
        for K in range(17,33):
            wl=self.common[0]+self.table[-1][K][0];wu=self.common[1]+self.table[-1][K][1]
            vl,vu=self.root(wl,wu);self.root_requests+=1
            prob=comb(32,K)*p**K*(self.S-p)**(32-K)
            pp=self.p[1][-1];probU=comb(32,K)*pp**K*(self.S-pp)**(32-K)
            tailMass[0]+=prob;tailMass[1]+=probU;term=prob*vu;tailU+=term
            tail.append({'K':K,'W_units':[str(wl),str(wu)],'v_units':[str(vl),str(vu)],
                         'binomial_probability_numerator_by_law':[str(prob),str(probU)],'upper_weighted_numerator':str(term)})
        req(self.root_requests==735487,'registered prefix and all16tail root requests')
        req(all(totalMass[i]+tailMass[i]==self.S**32 for i in range(2)),'both complete endpoint laws have exact unit mass')
        HL=Q(totalL,self.D);HU=Q(totalU+tailU,self.D);req(0<HL<=HU,'complete positive objective enclosure')
        return {'prefix_states':index,'prefix_count_by_K':countK,'state_bytes_sha256':globalhash.hexdigest(),
                'chunks':headers,'prefix_lower_at_pU_numerator':str(totalL),'prefix_upper_at_pL_numerator':str(totalU),
                'tail_upper_numerator':str(tailU),'tail_records':tail,'weighted_sum_denominator_power2':6336,
                'prefix_probability_numerator_by_law':list(map(str,totalMass)),
                'tail_probability_numerator_by_law':list(map(str,tailMass)),
                'probability_denominator_power2':6144,'both_endpoint_probability_laws_have_exact_unit_mass':True,
                'H_interval':[str(HL),str(HU)],'root_requests':self.root_requests,
                'root_cache_info':str(self.root.cache_info()),'root_cache_cap':65536,
                'tail_discarded':False,'prefix_renormalized':False,'CDF_endpoint_laws_independent_categories':False}
