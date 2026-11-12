from fractions import Fraction as F
for j in range(101):
    u = F(j,100)
    continuous = (2 - u**3 - (1-u)**3)/3
    assert continuous == F(1,3)+u-u*u
    assert ((1-u)**2+u**2)/4 + continuous/2 == F(5,12)
for m in range(2,33):
    for q in (1,2,3,5,9):
        exact = sum((1-F(abs(i-j),m))**q for i in range(m) for j in range(m))/m**2
        formula = F(1,m)+F(2,m)*sum((1-F(j,m))**(q+1) for j in range(1,m))
        assert formula == exact <= F(1,m)+F(2,q+2)
print('101 potential values and 155 exact posterior moment comparisons passed')
