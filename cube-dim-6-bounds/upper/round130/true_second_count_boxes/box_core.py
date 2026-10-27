"""Positive exact box bookkeeping, fixed Arb512 endpoint kernels."""
import gzip,hashlib,itertools,json,math,os
from fractions import Fraction as Q
from pathlib import Path
from functools import lru_cache

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(x,s):
    if not x:raise ArithmeticError(s)
def atomic(path,data):
    path=Path(path);tmp=path.with_name(path.name+'.own-tmp')
    with tmp.open('w')as f:json.dump(data,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)
def compressed(path,data):
    path=Path(path);tmp=path.with_name(path.name+'.own-tmp')
    payload=json.dumps(data,separators=(',',':'),allow_nan=False).encode()
    with tmp.open('wb')as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as f:f.write(payload)
        raw.flush();os.fsync(raw.fileno())
    os.replace(tmp,path)
def exact_ends(value):
    require(value.is_finite(),'finite Arb enclosure');out=[]
    for z in (value.lower(),value.upper()):
        m,e=z.man_exp();out.append(Q(int(m))*Q(2)**int(e))
    return out
def grid(value,bits):
    lo,hi=exact_ends(value);S=1<<bits
    return Q((lo*S).numerator//(lo*S).denominator,S),Q(-((-hi*S).numerator//(-hi*S).denominator),S)
def pair_sum(values):
    # Exact positive sum of all496 terms, no cancellation of interval sums.
    return sum((values[i]*values[j] for i in range(32)for j in range(i+1,32)),Q(0))
def read_cones(spec):
    import sympy as sp
    certificate=json.loads(Path(spec['audited_cells_path']).read_text());audit=json.loads(Path(spec['fan_audit_path']).read_text())
    R0=Q(spec['R0']);require(audit['status']=='PASS','audited fan PASS')
    require(Q(certificate['full_exp_minus_symbol_integral'])==Q(audit['full_exp_minus_symbol_integral'])==R0,'old audited R0 binding')
    cells=certificate['cells'];tokens=list(map(int,Path(spec['fan_path']).read_text().split()));require(len(cells)==tokens[0]==21,'21 source cones')
    V=[(1,)+tail for tail in itertools.product((-1,1),repeat=5)]
    require(sp.Matrix(V).T*sp.Matrix(V)==32*sp.eye(6),'tight original body32 frame')
    orbit={0:0,2:0,4:0}
    for i in range(32):
        for j in range(i+1,32):orbit[abs(sum(a*b for a,b in zip(V[i],V[j])))]+=1
    require(orbit=={4:96,2:240,0:160},'complete3 pair orbits')
    cursor=1;cones=[]
    for cid,cell in enumerate(cells):
        determinant=tokens[cursor];cursor+=1;rays=[];rates=[];scales=[]
        require(determinant>0,'positive source determinant')
        for j in range(6):
            rate=tokens[cursor];ray=tokens[cursor+1:cursor+7];cursor+=7
            audited=list(map(Q,cell['rays_a0_equals_one'][j]));scale=Q(ray[0])/audited[0]
            require(scale>0 and all(Q(a)==scale*b for a,b in zip(ray,audited)),'exact old ray binding')
            require(Q(rate)==scale*Q(cell['rates'][j]),'exact old rate binding')
            require(rate==sum(abs(sum(v[k]*ray[k]for k in range(6)))for v in V),'original psi rate')
            rays.append(ray);rates.append(rate);scales.append(scale)
        require(abs(int(sp.Matrix(rays).det()))==determinant,'source ray determinant')
        require(Q(determinant)==Q(cell['determinant'])*math.prod(scales),'audited source determinant scale')
        nodes=[]
        for v in V:
            dots=[sum(v[k]*r[k]for k in range(6))for r in rays]
            require(min(dots)>=0 or max(dots)<=0,'spectral sign consistency')
            b=[Q(abs(dot),rate)for dot,rate in zip(dots,rates)]
            require(any(b) and all(0<=z<Q(1,2)for z in b),'positive bounded actual spectral nodes')
            nodes.append([str(z)for z in b])
        cones.append({'id':cid,'alpha':str(Q(determinant,math.prod(rates))),'nodes_by_body':nodes,'rays':rays,'rates':rates,'determinant':determinant})
    require(cursor==len(tokens) and 46080*sum(Q(c['alpha'])for c in cones)==R0,'entire old fan normalization')
    return cones,orbit

class Evaluator:
    def __init__(self,spec):
        from flint import arb,fmpq,ctx
        ctx.prec=spec['precision_bits'];self.arb=arb;self.fmpq=fmpq;self.bits=spec['endpoint_grid_bits']
        self.h=Q(spec['h']);self.sh=self.A(self.h).sqrt();self.pL,self.pU=map(Q,spec['p_interval'])
        lo,hi=grid(self.sh.erf(),self.bits)
        require(self.pL<=lo<=hi<=self.pU,'fresh p enclosure inside pinned128 p interval')
        self.kernel_calls=0;self.vertex_calls=0;self.cone=None;self.b=None;self.vertex=self._vertex_cache()
    def A(self,q):q=Q(q);return self.arb(self.fmpq(q.numerator,q.denominator))
    def bind(self,cone):
        self.cone=cone;self.b=[[Q(z)for z in row]for row in cone['nodes_by_body']];self.vertex.cache_clear()
    def _vertex_cache(self):
        @lru_cache(maxsize=2048)
        def vertex(t):
            self.vertex_calls+=1;E=[-(1-self.A(z)).log()for z in t];out=[]
            for b in self.b:
                self.kernel_calls+=1
                if all(z==0 or t[j]==0 for j,z in enumerate(b)):
                    value=self.sh.erf()
                else:
                    a=sum((self.A(z)*E[j]for j,z in enumerate(b)),self.arb(0));u=a/(2*self.sh)
                    value=((u-self.sh).erfc()-(2*a).exp()*(u+self.sh).erfc())/2
                lo,hi=grid(value,self.bits)
                # Proven probabilistic bound, intersect only after full signed formula enclosure.
                lo=max(Q(0),lo);hi=min(self.pU,hi)
                require(0<=lo<=hi<=self.pU,'nonempty intersection with proven kernel range')
                out.append([str(lo),str(hi)])
            return out
        return vertex
    def box(self,lo,hi,prefix):
        lo=tuple(lo);hi=tuple(hi);require(len(lo)==len(hi)==6 and all(0<=a<b<1 for a,b in zip(lo,hi)),'honest rational6D box')
        kl=self.vertex(lo);ku=self.vertex(hi)
        low=pair_sum([Q(z[0])for z in ku]);up=pair_sum([Q(z[1])for z in kl]);require(0<=low<=up,'positive monotone corner bound')
        vol=math.prod(b-a for a,b in zip(lo,hi))
        return {'prefix':prefix,'lower_t':[str(z)for z in lo],'upper_t':[str(z)for z in hi],
                'volume':str(vol),'kernel_at_lower_t':kl,'kernel_at_upper_t':ku,
                'F_interval':[str(low),str(up)],'integral_interval':[str(vol*low),str(vol*up)]}
