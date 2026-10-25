% Public-model four-case MATLAB cross-check. Generator Delta.
pubU08a = [.5,.5,0.47716536955561023,0.47095317569629574];
pubU08n = [3,4,3,4];
for pubU08i=1:4
    pubU08dim=pubU08n(pubU08i);
    pubU08x=pubU08a(pubU08i);
    pubU08mass=integral(@(pubU08t) arrayfun(@(pubU08z) pubU08survival(pubU08z)^pubU08dim,pubU08t),0,0.12,'AbsTol',1e-12,'RelTol',1e-11) ...
        +integral(@(pubU08t) arrayfun(@(pubU08z) pubU08survival(pubU08z)^pubU08dim,pubU08t),0.12,Inf,'AbsTol',1e-12,'RelTol',1e-11);
    pubU08green=integral(@(pubU08t) arrayfun(@(pubU08z) pubU08heat(pubU08z,pubU08x)^pubU08dim,pubU08t),0,0.12,'AbsTol',1e-12,'RelTol',1e-11) ...
        +integral(@(pubU08t) arrayfun(@(pubU08z) pubU08heat(pubU08z,pubU08x)^pubU08dim,pubU08t),0.12,Inf,'AbsTol',1e-12,'RelTol',1e-11);
    pubU08cost=pubU08mass/(2^pubU08dim*pubU08x^pubU08dim*pubU08green);
    fprintf('PUBLIC_GREEN n=%d a=%.17g R=%.17g mass=%.17g green=%.17g cost=%.17g\n',pubU08dim,pubU08x,1/pubU08x,pubU08mass,pubU08green,pubU08cost);
end

function pubU08value=pubU08heat(pubU08t,pubU08x)
if pubU08t==0
    pubU08value=0;
elseif pubU08t<0.12
    pubU08m=-4:4;
    pubU08value=sum(exp(-(pubU08x-4*pubU08m).^2/(4*pubU08t))-exp(-(pubU08x-(4*pubU08m+2)).^2/(4*pubU08t)))/sqrt(4*pi*pubU08t);
else
    pubU08odd=1:2:79;
    pubU08value=sum(exp(-(pi*pubU08odd/2).^2*pubU08t).*cos(pi*pubU08odd*pubU08x/2));
end
end

function pubU08value=pubU08survival(pubU08t)
if pubU08t==0
    pubU08value=1;
elseif pubU08t<0.12
    pubU08k=0:4;
    pubU08value=1-2*sum((-1).^pubU08k.*erfc((2*pubU08k+1)/(2*sqrt(pubU08t))));
else
    pubU08k=0:39;
    pubU08odd=2*pubU08k+1;
    pubU08value=4/pi*sum((-1).^pubU08k./pubU08odd.*exp(-(pi*pubU08odd/2).^2*pubU08t));
end
end
