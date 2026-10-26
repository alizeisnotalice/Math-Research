"""111 pure mathematics core. No data generation at import time."""
from fractions import Fraction as Q
import hashlib,itertools

def poly(p,u):return p[0]+p[1]*u+p[2]*u*u
def quad_range(p,a,b):
 pts=[a,b]
 if p[2]:
  v=-p[1]/(2*p[2])
  if a<v<b:pts.append(v)
 vals=[(poly(p,u),u)for u in pts];lo=min(x[0]for x in vals);hi=max(x[0]for x in vals)
 return lo,hi,[u for z,u in vals if z==lo],[u for z,u in vals if z==hi]
def interval(lo,hi):return {'lower':str(lo),'upper':str(hi)}

class Arithmetic:
 def __init__(self):
  from flint import arb,fmpq,ctx
  ctx.prec=384;self.arb=arb;self.fmpq=fmpq
 def A(self,x):
  x=Q(x);return self.arb(self.fmpq(x.numerator,x.denominator))
 def end(self,x,side):
  m,e=(x.lower()if side=='lower'else x.upper()).man_exp();return Q(int(m))*Q(2)**int(e)
 def enc(self,x):
  l,e=x.lower().man_exp();u,f=x.upper().man_exp();e=int(e);f=int(f);g=min(e,f)
  return {'lower_mantissa':str(int(l)*2**(e-g)),'upper_mantissa':str(int(u)*2**(f-g)),'common_binary_exponent':str(g)}

def generate_geometry(k):
 """Called once per new k by parent after release, never by workers."""
 assert k in(2,3,4,5);m=6-k;ar=Arithmetic();A=ar.A
 V=[(1,)+tuple(z)for z in itertools.product([-1,1],repeat=5)]
 excluded={'v':31,'w':V.index(tuple([1]*k+[-1]*m))};branches={}
 for label,omit in excluded.items():
  table={};digest=hashlib.sha256();visits=0
  for i in range(32):
   if i==omit:continue
   for j in range(32):
    if j==omit or i==j:continue
    ai=sum(V[i][:k]);bi=sum(V[i][k:]);aj=sum(V[j][:k]);bj=sum(V[j][k:])
    c=sum(x*y for x,y in zip(V[i],V[j]));D=36-c*c;assert D in(20,32,36)
    p=(Q(k)-Q(6*ai*ai-2*c*ai*aj+6*aj*aj,D),-Q(12*ai*bi-2*c*(ai*bj+aj*bi)+12*aj*bj,D),Q(m)-Q(6*bi*bi-2*c*bi*bj+6*bj*bj,D))
    table[(p,D)]=table.get((p,D),0)+1;visits+=1
    digest.update((f'{i},{j};{p[0]},{p[1]},{p[2]};{D}\n').encode())
  assert visits==930 and sum(table.values())==930
  hist={D:sum(n for(p,d),n in table.items()if d==D)for D in(20,32,36)}
  assert hist=={20:180,32:450,36:300}
  assert len(table)<=((k*(m+1))*(k*(m+1)+1)//2)*5
  rows=[]
  for idx,((p,D),n)in enumerate(sorted(table.items(),key=lambda x:(x[0][0],x[0][1]))):
   lo,hi,atlo,athi=quad_range(p,Q(1,4),Q(3,4));assert lo>Q(k,4096)
   gap=A(lo).sqrt()-A(k).sqrt()/64;assert gap.lower()>0
   rows.append({'class':f'k{k}{label}{idx:03d}','P':[str(x)for x in p],'J2_squared':D,'n':n,'P_min_exact':str(lo),'P_max_exact':str(hi),'min_nodes':[str(x)for x in atlo],'max_nodes':[str(x)for x in athi],'full_radius_positive_gap':ar.enc(gap)})
  branches[label]={'excluded_unoriented_index':omit,'ordered_visits':930,'class_count':len(rows),'count_by_J2_squared':hist,'pair_visit_digest':digest.hexdigest(),'classes':rows}
 return {'k':k,'m':m,'u':['1/4','3/4'],'radius':'1/64','aggregation':'equal exact(P,J2) only; not geometric orbit claim','branches':branches}

def clamp_breaks(k,m,a,b):
 pts={a,b}
 for t in(Q(-1,2),Q(1,2)):
  x=(6*t-k)/m
  if a<x<b:pts.add(x)
 for t in(Q(-1,16),Q(1,16)):
  x=(k-6*t)/m
  if a<x<b:pts.add(x)
 if a<Q(1,2)<b:pts.add(Q(1,2))
 return sorted(pts)

class CellEngine:
 def __init__(self,job,geometry):
  self.job=job;self.k=int(job['k']);self.m=6-self.k;self.ar=Arithmetic();A=self.ar.A
  self.delta=Q(1,64);self.q=Q(6,5);self.sqrtk=A(self.k).sqrt();self.sqrtkm6=A(Q(self.k*self.m,6)).sqrt()
  self.sqrtD={D:A(D).sqrt()for D in(20,32,36)}
  ap=1/(2*A(self.k*self.m).sqrt());self.aplo=self.ar.end(ap,'lower');self.aphi=self.ar.end(ap,'upper')
  pi=self.ar.arb.pi();v4=pi*pi/2
  mass=A(Q(10,7))*v4*A(Q(self.m,self.k)).sqrt()*A(self.delta)**A(Q(14,5))
  self.masslo=self.ar.end(mass,'lower');self.masshi=self.ar.end(mass,'upper')
  self.groups={b:[(x['class'],tuple(Q(z)for z in x['P']),int(x['J2_squared']),int(x['n']))for x in geometry['branches'][b]['classes']]for b in('v','w')}

 def origin_poly(self,label,a,b):
  k,m=self.k,self.m;t=Q(1,2)if label=='v'else Q(1,16)
  star=lambda u:(k+m*u)/6 if label=='v'else(k-m*u)/6
  lo=min(star(a),star(b));hi=max(star(a),star(b))
  if lo>=t:fixed=t
  elif hi<=-t:fixed=-t
  else:
   assert lo>=-t and hi<=t,'missing actual clamp breakpoint'
   c=Q(k*m,6);return(c,-2*c if label=='v' else 2*c,c)
  return(Q(k)*(1-fixed)**2+m*fixed**2,-2*m*fixed if label=='v' else 2*m*fixed,Q(m))

 def branch_range(self,label,a,b,rlo,rhi):
  A=self.ar.A;Wlo=A(0);Whi=A(0);gaps=[]
  for name,p,D,n in self.groups[label]:
   pmin,pmax,_,_=quad_range(p,a,b)
   denlo=A(pmin).sqrt()-self.sqrtk*A(rhi);denhi=A(pmax).sqrt()-self.sqrtk*A(rlo)
   assert denlo.lower()>0 and denhi.lower()>0
   gaps.append(self.ar.end(denlo,'lower'))
   Whi+=A(n)/(self.sqrtD[D]*denlo);Wlo+=A(n)/(self.sqrtD[D]*denhi)
  wl=self.ar.end(Wlo,'lower');wu=self.ar.end(Whi,'upper')
  fmin,fmax,_,_=quad_range(self.origin_poly(label,a,b),a,b)
  h0l=self.ar.end(A(fmin).sqrt()-self.sqrtk*A(rhi),'lower');h0u=self.ar.end(A(fmax).sqrt()-self.sqrtk*A(rlo),'upper')
  if label=='v':
   if b<=Q(1,2):ml,mu=a,b
   else:assert a>=Q(1,2);ml,mu=1-b,1-a
  else:ml,mu=Q(7,8)-b,Q(7,8)-a
  h1l=self.ar.end(self.sqrtkm6*A(ml)-self.sqrtk*A(rhi),'lower');h1u=self.ar.end(self.sqrtkm6*A(mu)-self.sqrtk*A(rlo),'upper')
  assert min(h0l,h1l)>0
  lower=self.ar.end(A(wl)/A(h0u*h0u*h1u*h1u),'lower');upper=self.ar.end(A(wu)/A(h0l*h0l*h1l*h1l),'upper')
  return lower,upper,min(gaps+[h0l,h1l])

 def cell(self,state,a,b,s0,s1,du,ds,cid,seed):
  A=self.ar.A;state['cell_bounds']+=1
  rlo=self.delta*s0**5;rhi=self.delta*s1**5
  vl,vu,vg=self.branch_range('v',a,b,rlo,rhi);wl,wu,wg=self.branch_range('w',a,b,rlo,rhi)
  tl=self.aplo+rlo*(vl+wl);tu=self.aphi+rhi*(vu+wu)
  pl=self.ar.end(A(tl)**A(self.q),'lower');pu=self.ar.end(A(tu)**A(self.q),'upper')
  assert pl>0 and pl**5<=tl**6 and pu**5>=tu**6
  weight=(b-a)*(s1**14-s0**14);L=self.masslo*weight*pl;U=self.masshi*weight*pu;assert 0<=L<=U
  return {'id':cid,'seed':seed,'a':a,'b':b,'s0':s0,'s1':s1,'u_depth':du,'s_depth':ds,'lower':L,'upper':U,'T_lower':tl,'T_upper':tu,'power_lower':pl,'power_upper':pu,'weight':weight,'positive_gap_lower':min(vg,wg)}

 def initial_state(self):
  import heapq
  a,b=map(Q,self.job['u']);ub=clamp_breaks(self.k,self.m,a,b);sb=[Q(0),Q(1,2),Q(1)]
  st={'job':self.job,'initial_u_partition':ub,'initial_seeds':[],'refinements':[],'leaves':{},'heap':[],'nextid':0,'cell_bounds':0,'splits':0,'lower':Q(0),'upper':Q(0),'initial_leaves':0,'mass_constant':interval(self.masslo,self.masshi),'large_shape_coefficient_1_over_Jvw':interval(self.aplo,self.aphi),'completed_phases':[]}
  for ua,uz in zip(ub,ub[1:]):
   for sa,sz in zip(sb,sb[1:]):
    cid=st['nextid'];st['nextid']+=1;c=self.cell(st,ua,uz,sa,sz,0,0,cid,cid)
    st['leaves'][cid]=c;heapq.heappush(st['heap'],(-(c['upper']-c['lower']),cid));st['initial_seeds'].append({'id':cid,'u':[ua,uz],'s':[sa,sz]});st['lower']+=c['lower'];st['upper']+=c['upper']
  st['initial_leaves']=len(st['leaves']);return st

 def run_to(self,st,limit,depth,stop_requested,progress):
  import heapq
  before={'leaves':len(st['leaves']),'cell_bounds':st['cell_bounds'],'refinements':len(st['refinements']),'nextid':st['nextid']};outcome=None
  st['phase_entry_counts']=before
  maxbounds=2*limit-st['initial_leaves']
  while True:
   if stop_requested():outcome='STOP_SAVE_REQUESTED';break
   if len(st['leaves'])>=limit:outcome='FINITE_LEAF_CAP';break
   if st['cell_bounds']+2>maxbounds:outcome='FINITE_CELL_BOUND_CAP';break
   chosen=None
   while st['heap']:
    _,cid=heapq.heappop(st['heap'])
    if cid not in st['leaves']:continue
    c=st['leaves'][cid]
    if c['u_depth']+c['s_depth']>=depth:continue
    chosen=c;break
   if chosen is None:outcome='FINITE_DEPTH_CAP';break
   a,b,s0,s1=chosen['a'],chosen['b'],chosen['s0'],chosen['s1'];du,ds=chosen['u_depth'],chosen['s_depth'];axis='u' if du<=ds else 's'
   ids=[st['nextid'],st['nextid']+1];st['nextid']+=2
   if axis=='u':
    mid=(a+b)/2;children=[self.cell(st,a,mid,s0,s1,du+1,ds,ids[0],chosen['seed']),self.cell(st,mid,b,s0,s1,du+1,ds,ids[1],chosen['seed'])]
   else:
    mid=(s0+s1)/2;children=[self.cell(st,a,b,s0,mid,du,ds+1,ids[0],chosen['seed']),self.cell(st,a,b,mid,s1,du,ds+1,ids[1],chosen['seed'])]
   del st['leaves'][chosen['id']]
   for c in children:st['leaves'][c['id']]=c;heapq.heappush(st['heap'],(-(c['upper']-c['lower']),c['id']))
   st['lower']+=sum((x['lower']for x in children),Q(0))-chosen['lower'];st['upper']+=sum((x['upper']for x in children),Q(0))-chosen['upper']
   st['refinements'].append({'parent':chosen['id'],'axis':axis,'children':ids});st['splits']+=1
   if st['splits']%256==0:progress(st)
  assert st['lower']==sum((x['lower']for x in st['leaves'].values()),Q(0))and st['upper']==sum((x['upper']for x in st['leaves'].values()),Q(0))
  assert len(st['leaves'])<=limit and st['cell_bounds']<=maxbounds
  st['last_outcome']=outcome;st['phase_entry_counts']=before
  return st

 def snapshot(self,st):
  return{k:v for k,v in st.items()if k!='heap'}
