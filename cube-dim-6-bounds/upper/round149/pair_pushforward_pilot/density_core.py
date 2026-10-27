"""Exact rational rank2 compiler for a finite released pattern; no k evaluation."""
from fractions import Fraction as F
from math import factorial

def need(v,msg):
 if not v:raise ArithmeticError(msg)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def det(a,b):return a[0]*b[1]-a[1]*b[0]
def direction(v):return tuple(z/sum(v)for z in v)

def compile_case(case,caps,checkpoint,counter_state):
 import sympy as sp
 need(sp.__version__=='1.14.0','PINNED_SYMPY')
 X,Y,R,U,V,S,T=sp.symbols('X Y R U V S T')
 counts=counter_state;counts.update({'density_stage_check_attempts':0,'density_stage_checks':0,'convolution_term_branches':0,'resonant_branches':0,'nonresonant_branches':0,'literal_endpoint_substitutions':0,'Laplace_monomial_summands':0,'checked_expression_node_visits':0,'expression_nodes_max':0,'nonzero_terms_max':0,'polynomial_monomials_max':0,'polynomial_degree_max':0,'coefficient_bits_max':0,'new_k_integral_special_function_native_kernel_calls':0,'old_6D_grid_replays':0})
 def q(z):return sp.Rational(z.numerator,z.denominator)
 def linear(z):return q(z[0])*X+q(z[1])*Y
 def tree(expr):
  n=0;stack=[expr]
  while stack:
   need(n<caps['expression_nodes_per_checked_expression'],'EXPRESSION_NODE_CAP');need(counts['checked_expression_node_visits']<caps['checked_expression_node_visits_per_case'],'CHECKED_NODE_VISIT_CAP')
   z=stack.pop();n+=1;counts['checked_expression_node_visits']+=1;counts['expression_nodes_max']=max(counts['expression_nodes_max'],n);stack.extend(z.args)
 def rational_bits(z):
  z=sp.Rational(z);n=max(abs(int(z.p)).bit_length(),int(z.q).bit_length());counts['coefficient_bits_max']=max(counts['coefficient_bits_max'],n);need(n<=caps['rational_coefficient_bits'],'RATIONAL_COEFFICIENT_BIT_CAP')
 def poly(expr):
  p=sp.Poly(sp.expand(expr),X,Y,domain=sp.QQ)
  for c in p.coeffs():rational_bits(c)
  return p.as_expr()
 def bounds(left,right,v,sample):
  lower=[(F(0),F(0))];upper=[]
  for L in [(-left[1],left[0]),(right[1],-right[0])]:
   d=dot(L,v)
   if not d:
    if dot(L,sample)<0:return None
   elif d>0:upper.append(tuple(z/d for z in L))
   else:lower.append(tuple(z/d for z in L))
  need(bool(upper),'POINTED_CONE_FINITE_UPPER')
  a=max(lower,key=lambda z:dot(z,sample));b=min(upper,key=lambda z:dot(z,sample))
  return None if dot(b,sample)<=dot(a,sample)else(a,b)
 def verify(cells,active):
  need(counts['density_stage_check_attempts']<caps['density_stage_checks_per_case'],'STAGE_CHECK_CAP');counts['density_stage_check_attempts']+=1
  dirs={direction(v)for v in active};need(len(cells)<=len(dirs)-1<=5,'DISTINCT_DIRECTION_FAN_BOUND')
  nterms=sum(len(terms)for _,_,terms in cells);counts['nonzero_terms_max']=max(counts['nonzero_terms_max'],nterms);need(nterms<=caps['nonzero_terms_per_stage'],'DENSITY_TERM_CAP')
  summands=[];encoded=[];monos=0
  for a,b,terms in cells:
   D=det(a,b);need(D>0,'ORIENTED_SECTOR');rows=[]
   for lam,p in sorted(terms.items()):
    pp=sp.Poly(p,X,Y,domain=sp.QQ);degree=pp.total_degree();counts['polynomial_degree_max']=max(counts['polynomial_degree_max'],degree);need(degree<=len(active)-2<=4,'POLYNOMIAL_DEGREE_BOUND')
    monos+=len(pp.terms());need(monos<=caps['polynomial_monomials_per_stage'],'POLYNOMIAL_MONOMIAL_CAP');tree(p)
    for z in lam:rational_bits(q(z))
    xy=sp.Poly(sp.expand(p.subs({X:q(a[0])*U+q(b[0])*V,Y:q(a[1])*U+q(b[1])*V},simultaneous=True)),U,V,domain=sp.QQ)
    L=(q(lam[0])+S)*q(a[0])+(q(lam[1])+T)*q(a[1]);W=(q(lam[0])+S)*q(b[0])+(q(lam[1])+T)*q(b[1])
    for (i,j),coef in xy.terms():
     summands.append(q(D)*coef*factorial(i)*factorial(j)/(L**(i+1)*W**(j+1)));need(len(summands)<=caps['Laplace_monomial_summands_per_stage'],'LAPLACE_SUMMAND_CAP')
    rows.append({'lambda':list(map(str,lam)),'polynomial':str(p)})
   encoded.append({'left_ray':list(map(str,a)),'right_ray':list(map(str,b)),'terms':rows})
  counts['polynomial_monomials_max']=max(counts['polynomial_monomials_max'],monos);counts['Laplace_monomial_summands']+=len(summands)
  total=sp.Add(*summands,evaluate=False);tree(total);normal=sp.cancel(total);tree(normal)
  expected=sp.prod(1/(1+q(v[0])*S+q(v[1])*T)for v in active)
  difference=sp.cancel(normal-expected);need(difference==0,'EXACT_LAPLACE_DIFFERENCE_NUMERATOR_ZERO')
  need(sp.simplify(normal.subs({S:0,T:0}))==1,'WHOLE_NORMALIZATION_ONE');counts['density_stage_checks']+=1
  stage={'active_columns':len(active),'active_column_vectors':[list(map(str,z))for z in active],'distinct_directions':len(dirs),'open_sector_count':len(cells),'nonzero_terms':nterms,'polynomial_monomials':monos,'Laplace_monomial_summands':len(summands),'density_cells':encoded,'Laplace_rational_normal_form':str(normal),'Laplace_expected_product':str(expected),'Laplace_difference_numerator':'0','normalization':1,'positivity':'complete density is a convolution of probability laws; signed terms are not separately positive','actual_counts_to_this_stage':dict(counts)}
  checkpoint(stage)
 columns=[tuple(map(F,v))for v in case['columns']];need(len(columns)==6 and all(min(v)>=0 for v in columns),'SIX_NONNEGATIVE_COLUMNS')
 cols=[v for v in columns if any(v)];m=len(cols);need(2<=m<=6,'ACTIVE_M2_TO_M6')
 pivot=next((i,j)for i in range(m)for j in range(i+1,m)if det(cols[i],cols[j]))
 u,w=cols[pivot[0]],cols[pivot[1]]
 if det(u,w)<0:u,w=w,u
 D=det(u,w);lam=((w[1]-u[1])/D,(u[0]-w[0])/D);need(dot(lam,u)==dot(lam,w)==1,'INITIAL_COVECTOR')
 a,b=direction(u),direction(w);cells=[(a,b,{lam:q(1/D)})];fan=[a,b];active=[u,w];verify(cells,active)
 for index,v in enumerate(cols):
  if index in pivot:continue
  newfan=sorted(set(fan+[direction(v)]),key=lambda z:z[1]);newcells=[]
  for a,b in zip(newfan,newfan[1:]):
   sample=tuple(a[j]+b[j]for j in range(2));newterms={}
   for left,right,terms in cells:
    interval=bounds(left,right,v,sample)
    if interval is None:continue
    lower,upper=interval
    for lam,p in terms.items():
     need(counts['literal_endpoint_substitutions']+2<=caps['literal_endpoint_substitutions_per_case'],'LITERAL_ENDPOINT_SUBSTITUTION_CAP');counts['convolution_term_branches']+=1
     eta=dot(lam,v)-1;pr=sp.Poly(sp.expand(p.subs({X:X-q(v[0])*R,Y:Y-q(v[1])*R},simultaneous=True)),R)
     if not eta:
      counts['resonant_branches']+=1;primitive=sum(c*R**(k[0]+1)/sp.Integer(k[0]+1)for k,c in pr.terms())
      # Count each literal substitution immediately before that call, even when the difference is merged once.
      counts['literal_endpoint_substitutions']+=1;high=primitive.subs(R,linear(upper));counts['literal_endpoint_substitutions']+=1;low=primitive.subs(R,linear(lower));value=poly(high-low);newterms[lam]=poly(newterms.get(lam,0)+value)
     else:
      counts['nonresonant_branches']+=1;expr=pr.as_expr();Q=sum((-1)**k*sp.diff(expr,R,k)/q(eta)**(k+1)for k in range(pr.degree()+1))
      for endpoint,sign in [(upper,1),(lower,-1)]:
       newlam=tuple(lam[j]-eta*endpoint[j]for j in range(2));counts['literal_endpoint_substitutions']+=1;value=poly(sign*Q.subs(R,linear(endpoint)));newterms[newlam]=poly(newterms.get(newlam,0)+value)
   newterms={lam:p for lam,p in newterms.items()if p!=0}
   if newterms:newcells.append((a,b,newterms))
  active.append(v);fan=newfan;cells=newcells;verify(cells,active)
 need(counts['density_stage_checks']==m-1,'EXACT_M_MINUS_ONE_CHECKS')
 need(counts['literal_endpoint_substitutions']==2*counts['convolution_term_branches'] and counts['convolution_term_branches']==counts['resonant_branches']+counts['nonresonant_branches'],'COUNTER_IDENTITIES')
 return {'status':'PASS_EXACT_DENSITY_CASE','case_id':case['case_id'],'canonical_key':case['canonical_key'],'active_columns':m,'distinct_directions':len({direction(v)for v in cols}),'actual_counts':counts,'new_k_integral_evaluations':0,'new_beta2_upper':False,'new_C6_bound':False}
