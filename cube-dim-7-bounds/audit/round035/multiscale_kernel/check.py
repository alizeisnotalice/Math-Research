from fractions import Fraction as Q
from sympy import symbols, integrate, simplify
w,h,x=symbols('w h x', positive=True)
# For 0<h<w/2, symmetry identifies a single boundary layer of width h at each end.
area=simplify(2*integrate(h+x,(x,0,h))+integrate(2*h,(x,h,w-h)))
assert simplify(area-(2*h*w-h*h))==0
# At matching centered intervals of width h/2, every pair lies inside |u+v-d|<h.
# The joint product-measure probability and maximum row-degree are identically 1.
assert Q(1,2)+Q(1,2)<2
print('strip area:',area,'; normalized:',simplify(area/w**2))
print('two concentrated uniform marginals: probability=1, maximum degree=1')
rho,d=symbols('rho d', positive=True)
y=symbols('y', nonnegative=True)
assert simplify(integrate(y/rho,(y,0,d))-d**2/(2*rho))==0
assert simplify(integrate(y**2/rho,(y,0,d))-d**3/(3*rho))==0
print('uniform-shell excess moments:',d**2/(2*rho),d**3/(3*rho))
