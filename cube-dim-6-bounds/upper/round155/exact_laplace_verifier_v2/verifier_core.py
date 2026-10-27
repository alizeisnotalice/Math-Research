"""Prepared independent saved-stage verifier; direct exact integer numerator, no SymPy."""
from fractions import Fraction as F
from math import gcd,lcm,factorial,comb
import ast,time
from api_ledger import Ledger

def need(v,m):
 if not v:raise ArithmeticError(m)
def dot(a,b):
 need(len(a)==len(b)==2,'DOT_DIMENSION_EXACTLY2');return sum((x*y for x,y in zip(a,b)),F(0))
def det(a,b):
 need(len(a)==len(b)==2,'DET_DIMENSION_EXACTLY2');return a[0]*b[1]-a[1]*b[0]
def direction(v):return tuple(x/sum(v)for x in v)

def verify_stage(stage,original_columns,caps,global_counts):
 # Import is intentionally inside the released Terra-only entry.
 import flint
 need(flint.__version__=='0.9.0','PINNED_PYTHON_FLINT')
 started=time.monotonic();cpu=time.process_time();k=stage['active_columns'];X=(1,0);Y=(0,1)
 ledger=Ledger(global_counts,caps,{'active_columns':k});call=ledger.invoke
 count={'rational_Gamma_terms':0,'native_terms_max':0,'native_coefficient_bits_max':0}
 def bitcheck(q):need(max(abs(q.numerator).bit_length(),q.denominator.bit_length())<=caps['rational_input_bits'],'RATIONAL_INPUT_BITS')
 def clean(p):
  p={v:c for v,c in p.items()if c};need(len(p)<=15 and all(sum(v)<=4 for v in p),'POLYNOMIAL_DEGREE_OR_MONOMIAL_CAP')
  for c in p.values():bitcheck(c)
  return p
 def add(p,q):
  r=dict(p)
  for v,c in q.items():r[v]=r.get(v,F(0))+c
  return clean(r)
 def mul(p,q):
  r={}
  for a,c in p.items():
   for b,d in q.items():
    v=(a[0]+b[0],a[1]+b[1]);r[v]=r.get(v,F(0))+c*d
  return clean(r)
 def parse(text):
  need(len(text)<=caps['polynomial_text_bytes'],'POLYNOMIAL_TEXT_CAP');tree=ast.parse(text,mode='eval');need(sum(1 for _ in ast.walk(tree))<=caps['polynomial_AST_nodes'],'POLYNOMIAL_AST_CAP')
  def f(n):
   if isinstance(n,ast.Constant)and type(n.value)is int:return clean({(0,0):F(n.value)})
   if isinstance(n,ast.Name)and n.id in ['X','Y']:return {X if n.id=='X'else Y:F(1)}
   if isinstance(n,ast.UnaryOp)and isinstance(n.op,(ast.USub,ast.UAdd)):
    p=f(n.operand);return clean({v:(-c if isinstance(n.op,ast.USub)else c)for v,c in p.items()})
   if isinstance(n,ast.BinOp):
    a=f(n.left);b=f(n.right)
    if isinstance(n.op,ast.Add):return add(a,b)
    if isinstance(n.op,ast.Sub):return add(a,{v:-c for v,c in b.items()})
    if isinstance(n.op,ast.Mult):return mul(a,b)
    if isinstance(n.op,ast.Div):
     need(set(b)=={(0,0)}and b[(0,0)]!=0,'DIVISION_ONLY_BY_NONZERO_RATIONAL');return clean({v:c/b[(0,0)]for v,c in a.items()})
    if isinstance(n.op,ast.Pow):
     need(set(b)<={(0,0)},'INTEGER_POWER');z=b.get((0,0),F(0));need(z.denominator==1 and 0<=z<=4,'POWER_CAP');r={(0,0):F(1)}
     for _ in range(int(z)):r=mul(r,a)
     return r
   raise ArithmeticError('UNSUPPORTED_POLYNOMIAL_AST')
  return f(tree.body)
 def expand_rays(p,a,b):
  out={}
  for (i,j),c in p.items():
   for h in range(i+1):
    for z in range(j+1):
     exponent=(h+z,i+j-h-z);value=c*comb(i,h)*comb(j,z)*a[0]**h*b[0]**(i-h)*a[1]**z*b[1]**(j-z)
     out[exponent]=out.get(exponent,F(0))+value
  return clean(out)
 def primitive(affine):
  # affine=(S coefficient,T coefficient,constant), primitive first slope positive.
  for v in affine:bitcheck(v)
  den=lcm(*(v.denominator for v in affine));ints=tuple(int(v*den)for v in affine);need(max(map(lambda z:abs(z).bit_length(),ints))<=caps['rational_input_bits'],'PRIMITIVE_FACTOR_BITS');g=gcd(*ints);need(g>0,'NONZERO_AFFINE_FACTOR')
  prim=tuple(v//g for v in ints);first=next(v for v in prim if v)
  if first<0:prim=tuple(-v for v in prim);g=-g
  return prim,F(g,den)
 terms=[];raw_affines=[]
 def append(coefficient,denominators):
  powers={}
  for affine,e in denominators:
   raw_affines.append(affine);factor,scale=primitive(affine);coefficient/=scale**e
   if factor==(0,0,1):continue
   powers[factor]=powers.get(factor,0)+e
  bitcheck(coefficient)
  if coefficient:terms.append((coefficient,powers))
  need(len(terms)<=caps['rational_terms_per_stage'],'RATIONAL_TERM_CAP')
 need(type(k)is int,'INTEGER_ACTIVE_COLUMN_COUNT')
 need(isinstance(original_columns,(list,tuple))and len(original_columns)==6 and all(isinstance(v,(list,tuple))and len(v)==2 for v in original_columns),'EXACT_RAW_COLUMN_DIMENSIONS')
 cols=[]
 for raw in original_columns:
  v=tuple(map(F,raw))
  for value in v:bitcheck(value)
  need(min(v)>=0,'NONNEGATIVE_RAW_COLUMNS');cols.append(v)
 cols=[v for v in cols if any(v)]
 pivot=next((i,j)for i in range(len(cols))for j in range(i+1,len(cols))if det(cols[i],cols[j]))
 ordered=[cols[i]for i in pivot]+[v for i,v in enumerate(cols)if i not in pivot];active=ordered[:k]
 need(2<=k<=len(cols)<=6 and stage['active_column_vectors']is not None,'ACTIVE_PREFIX')
 # Pivot orientation may swap its two vectors without altering their product law.
 need(len(stage['active_column_vectors'])==k and all(isinstance(v,(list,tuple))and len(v)==2 for v in stage['active_column_vectors']),'EXACT_SAVED_ACTIVE_DIMENSIONS')
 saved_active=[tuple(map(F,v))for v in stage['active_column_vectors']]
 for v in saved_active:
  for value in v:bitcheck(value)
 need(sorted(saved_active)==sorted(active),'UNCHANGED_SOURCE_ACTIVE_COLUMNS')
 rays=sorted({direction(v)for v in active},key=lambda v:v[1]);expected=list(zip(rays,rays[1:]));cells=stage['density_cells'];need(len(cells)==len(expected),'EXACT_SUPPORT_FAN')
 need(len(cells)<=5 and sum(len(c['terms'])for c in cells)<=2000,'SOURCE_DENSITY_SIZE')
 for cell,(a,b)in zip(cells,expected):
  need(isinstance(cell['left_ray'],(list,tuple))and isinstance(cell['right_ray'],(list,tuple))and len(cell['left_ray'])==len(cell['right_ray'])==2,'EXACT_RAY_DIMENSIONS')
  saved_a=tuple(map(F,cell['left_ray']));saved_b=tuple(map(F,cell['right_ray']))
  for value in saved_a+saved_b:bitcheck(value)
  need(saved_a==a and saved_b==b and det(a,b)>0,'ORIENTED_SOURCE_SUPPORT')
  for term in cell['terms']:
   need(isinstance(term['lambda'],(list,tuple))and len(term['lambda'])==2,'EXACT_LAMBDA_DIMENSION')
   lam=tuple(map(F,term['lambda']))
   for value in lam:bitcheck(value)
   p=parse(term['polynomial']);need(p and max(map(sum,p))<=k-2,'SOURCE_POLYNOMIAL_DEGREE')
   xy=expand_rays(p,a,b)
   for (i,j),c in xy.items():
    append(det(a,b)*c*factorial(i)*factorial(j),[((a[0],a[1],dot(lam,a)),i+1),((b[0],b[1],dot(lam,b)),j+1)])
 target=[((v[0],v[1],F(1)),1)for v in active];append(F(-1),target)
 # All nonconstant raw denominator slopes are nonnegative with positive sum.
 tilt=0
 for a,b,c in raw_affines:
  need(a>=0 and b>=0,'NONNEGATIVE_RAY_SLOPES')
  if a+b:tilt=max(tilt,max(0,int((-c)//(a+b))+1))
  else:need(c>0,'POSITIVE_CONSTANT_FACTOR')
 need(all(a*tilt+b*tilt+c>0 for a,b,c in raw_affines),'STRICT_POSITIVE_TILT')
 maximum={}
 for coefficient,powers in terms:
  for factor,e in powers.items():maximum[factor]=max(maximum.get(factor,0),e)
 need(len(maximum)<=caps['common_factor_count'],'COMMON_FACTOR_COUNT_CAP');degree=sum(maximum.values());need(degree<=caps['common_degree'],'COMMON_DENOMINATOR_DEGREE_CAP')
 B=1
 for coefficient,_ in terms:B=lcm(B,coefficient.denominator);need(B.bit_length()<=caps['scalar_LCM_bits'],'SCALAR_LCM_BIT_CAP')
 ctx=call('setup','fmpz_mpoly_ctx.get',lambda:flint.fmpz_mpoly_ctx.get(('S','T'),ordering='lex'))
 def inspect(p):
  length=call('query','fmpz_mpoly.__len__',lambda:len(p));zero=call('query','fmpz_mpoly.is_zero',p.is_zero)
  need(length<=caps['native_polynomial_terms'],'NATIVE_POLYNOMIAL_TERM_CAP')
  if not zero:need(call('query','fmpz_mpoly.total_degree',p.total_degree)<=degree,'NATIVE_POLYNOMIAL_DEGREE_CAP')
  count['native_terms_max']=max(count['native_terms_max'],length)
  coefficients=call('query','fmpz_mpoly.coeffs',p.coeffs);bits=0
  for coefficient in coefficients:
   value=call('query','fmpz.__int__',lambda coefficient=coefficient:int(coefficient));bits=max(bits,abs(value).bit_length())
  need(bits<=caps['native_coefficient_bits'],'NATIVE_COEFFICIENT_BITS');count['native_coefficient_bits_max']=max(count['native_coefficient_bits_max'],bits);return p
 def constant(z):return call('construction','fmpz_mpoly_ctx.constant',lambda:ctx.constant(z))
 factor_powers={}
 for f,e in maximum.items():
  data={(1,0):f[0],(0,1):f[1],(0,0):f[2]};base=call('construction','fmpz_mpoly_ctx.from_dict',lambda data=data:ctx.from_dict(data));inspect(base);factor_powers[f]=[constant(1)]
  for n in range(1,e+1):
   previous=factor_powers[f][-1];factor_powers[f].append(inspect(call('arithmetic','fmpz_mpoly.__mul__',lambda previous=previous,base=base:previous*base)))
 numerator=constant(0)
 for coefficient,powers in terms:
  scaled=coefficient*B;need(scaled.denominator==1,'INTEGER_SCALAR_CLEARING');product=constant(scaled.numerator)
  for factor,e in maximum.items():
   missing=e-powers.get(factor,0)
   if missing:
    right=factor_powers[factor][missing];product=inspect(call('arithmetic','fmpz_mpoly.__mul__',lambda product=product,right=right:product*right))
  numerator=inspect(call('arithmetic','fmpz_mpoly.__add__',lambda numerator=numerator,product=product:numerator+product))
 need(call('query','fmpz_mpoly.is_zero.final',numerator.is_zero),'EXACT_INTEGER_NUMERATOR_NONZERO')
 ledger.check_quiescent()
 count['rational_Gamma_terms']=len(terms)-1
 return {'status':'PASS_DIRECT_EXACT_INTEGER_LAPLACE_NUMERATOR_ZERO','active_columns':k,'sectors':len(cells),'common_denominator_factors':[{'primitive_STC':list(map(str,f)),'maximum_power':e}for f,e in sorted(maximum.items())],'common_denominator_total_degree':degree,'scalar_LCM_bit_length':B.bit_length(),'strict_positive_tilt_integer':tilt,'exact_integer_numerator_is_zero':True,'normalization':1,'normalization_basis':'positive-tilted finite signed-measure Laplace uniqueness identifies original probability law; no per-term positivity or unshifted pole evaluation','size_counts':count,'API_stage_ledger':ledger.stage,'API_budget_semantics':'original3M cap covers all setup+construction+arithmetic attempts; query attempts/completions/exceptions are separately and fully recorded','wall_seconds':time.monotonic()-started,'CPU_seconds':time.process_time()-cpu,'new_density_generations':0,'k_or_special_function_integrals':0,'legacy_native_integral_kernel_calls':0}
