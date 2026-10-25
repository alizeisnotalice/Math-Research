# Stage2 deterministic choices before inspecting the Gauss pair

Stage1 selected only m=3 by the registered larger-ceiling rule. The next positive weighted measure is d*u^6 dR, with moments mu_j=E[d^(j+1)u^6], j0..3. Its two Gauss nodes are the roots of t²-A*t+B where A=(mu0*mu3-mu1*mu2)/(mu0*mu2-mu1²), B=(mu1*mu3-mu2²)/(mu0*mu2-mu1²). Use the nearest multiples of1/10000, verifying the selection interval and 0<r1<r2<1. No node optimization or retry grid.

For h(t)=t^(-1/2), hprime(ri)=-h(ri)/(2ri), so every cubic Hermite coefficient is a rational linear combination of h(r1),h(r2). Enclose each h value using rational endpoints of denominator10^14 and verify their squares around1/ri. Form each coefficient lower endpoint using the signs of those rational linear combinations; then ell(d)=d*sum coefficient_j*d^j is globally below sqrt(d). Coefficients may be negative.

With U=sum coefficient_j*mu_j and V fromStage1, choose a=floor(10^5*U/V)/10^5 if U>0. Verify the actual corner value cL*(3a²U-2a³V), using014 c_lower. If U<=0 or the actual corner does not exceed014 current required g, stop without spatial production. If it does exceed, report the precise result and record the fixed-a/adaptive spatial checkpoint before running it.
