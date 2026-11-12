from decimal import Decimal as D, getcontext
getcontext().prec = 45
w = D(1)/3
exact = (1-w*w)*(1+w).ln() + w*w*w.ln() + w
# 2*g*log(1+w/g) increases for g>0; its limit at zero is zero.
steps = 10000
right = sum(2*(D(i)/steps)*(1+w/(D(i)/steps)).ln() for i in range(1,steps+1))/steps
left = right - 2*(1+w).ln()/steps
assert left < exact < right
mu = D('1.2784')  # rigorous lower bound from round 234
assert mu - D(2).ln() - right > D('0.118')
print('F1(1/3)', exact, 'grid bracket', left, right)
