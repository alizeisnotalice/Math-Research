import sympy as s
u,t,a=s.symbols('u t a',positive=True,real=True)
# Squared two-center convolution by differentiating unequal Cauchy convolution.
h,k=s.symbols('h k',positive=True)
base=s.pi*(h+k)/(h*k*(a*a+(h+k)**2))
r=s.simplify(s.diff(base,h,k)/(4*h*k));r=s.factor(r.subs({h:t,k:t}))
p=s.factor(t**4*r/(4*s.pi**4))
print('p_t axis distance a centers:',p)
G=s.integrate(p,(t,0,s.oo));print('G:',s.simplify(G));assert s.simplify(G-3/(32*s.pi**3*a*a))==0
