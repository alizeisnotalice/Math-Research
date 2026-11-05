"""Independent exact audit of author's rational probability-family parameters."""
from fractions import Fraction as F
from pathlib import Path
import json, platform

cases=[]
for n in (1,2,3,8,32,128):
    q=F(4)**(-n)/16
    alpha=1+q
    a=[F(0)]+[F(2)**(j-n-1)/alpha for j in range(1,n+1)]
    epsilon=1-alpha**(-2)
    masses=[4*(a[j]**2-a[j-1]**2) for j in range(1,n+1)]
    assert sum(masses)+epsilon==1
    assert epsilon>0 and all(m>0 for m in masses)
    assert epsilon/a[1]**2==(2+q)/16<1
    t=[F(0)]+[epsilon/(16*a[j]) for j in range(1,n+1)]
    h=t[n]/4
    atom_loss=F(0)
    l1_loss=F(0)
    for j in range(1,n+1):
        assert sum(masses[:j])==4*a[j]**2
        rho=a[j]+t[j]
        assert 4*a[j]**2+epsilon-4*rho**2==epsilon/2-epsilon**2/(64*a[j]**2)>0
        assert h<=t[j]/4
        assert a[j]+h+t[j]/2<rho
        assert a[j]+h<rho
        assert h+t[j]/2<rho
        atom_loss+=a[j]*t[j]
        l1_loss+=a[j]*t[j]/2
    assert atom_loss/epsilon==F(n,16)
    assert l1_loss/epsilon==F(n,32)
    assert a[n]+h<F(1,2)
    assert 2*h+a[n]<1
    cases.append({'n':n,'q':str(q),'epsilon':str(epsilon),
                  'epsilon_over_a1_squared':str(epsilon/a[1]**2),
                  'mass_sum':str(sum(masses)+epsilon),
                  'atom_loss_over_epsilon':str(atom_loss/epsilon),
                  'l1_loss_over_epsilon':str(l1_loss/epsilon),
                  'right_support':str(a[n]+h),'h':str(h)})
out={'python':platform.python_version(),'arithmetic':'fractions.Fraction exact',
     'all_assertions_passed':True,'cases':cases,
     'scope':'Exact finite-parameter supporting arithmetic only; all-N proof is analytic.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
