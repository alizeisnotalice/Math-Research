from fractions import Fraction as F
from itertools import product
# Under independent permutations in each of two blocks, average an optimizer.
# Four orbit types (private or shared per block); buy in decreasing incidence multiplicity.
for N,L,t in [(4,16,F(1,4)),(8,32,F(1,8)),(17,64,F(1,10))]:
    for b in (1,2):
        required=t*L**b
        cost=F(0)
        orbits=[]
        for flags in product((0,1),repeat=b):
            private=sum(flags)
            within=(L-1)**private
            orbits.append((N**private,within))
        for rate,available in sorted(orbits):
            take=min(required,available)
            cost+=rate*take
            required-=take
        assert required==0 and cost>=N**b*L**b*(t-F(b,L))
        if b==1:
            assert cost==1+N*(t*L-1)
        print(N,L,b,t,cost)
