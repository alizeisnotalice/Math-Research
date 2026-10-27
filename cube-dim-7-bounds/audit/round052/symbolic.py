from sympy import symbols,Rational as Q,factor,cancel
A,k,t,N=symbols('A k t N',positive=True)
a=A/8;d=k/N;s=(4*A+t)/16;P=Q(3579,4096);eps=Q(1,1024)
mass=Q(63,64)*d/(2*(1+3*d/4));br=(1-3/a)*P/(2*(1+eps))-2*eps
c=(1-eps)*mass*br/4-s*d**2/(1-s*d**2)
fnum=512*3579*(A-24)-8200*A;pn=(4*A+t)*k*k;pd=16*N*N-pn;n=64449*k*fnum;z=131072*(4*N+3*k)*4198400*A
assert cancel(c-(n*pd-pn*z)/(z*pd))==0
assert cancel((1-eps)*mass*br/4-n/z)==0
print('SymPy exact symbolic identity: residual',factor(c-(n*pd-pn*z)/(z*pd)))
print('first term:',factor((1-eps)*mass*br/4))
