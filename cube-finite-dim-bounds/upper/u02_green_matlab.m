% U02 candidate-only independent MATLAB check. Generator is Delta.
upU02a = [0.47716536955561023,0.47095317569629574];
upU02n = [3,4];
for upU02i=1:2
    upU02dim=upU02n(upU02i);
    upU02x=upU02a(upU02i);
    upU02mass=integral(@(upU02t) arrayfun(@(upU02z) upU02survival(upU02z)^upU02dim,upU02t),0,0.12,'AbsTol',1e-12,'RelTol',1e-11) ...
        +integral(@(upU02t) arrayfun(@(upU02z) upU02survival(upU02z)^upU02dim,upU02t),0.12,Inf,'AbsTol',1e-12,'RelTol',1e-11);
    upU02green=integral(@(upU02t) arrayfun(@(upU02z) upU02heat(upU02z,upU02x)^upU02dim,upU02t),0,0.12,'AbsTol',1e-12,'RelTol',1e-11) ...
        +integral(@(upU02t) arrayfun(@(upU02z) upU02heat(upU02z,upU02x)^upU02dim,upU02t),0.12,Inf,'AbsTol',1e-12,'RelTol',1e-11);
    upU02cost=upU02mass/(2^upU02dim*upU02x^upU02dim*upU02green);
    fprintf('U02 n=%d a=%.17g R=%.17g mass=%.17g green=%.17g cost=%.17g\n',upU02dim,upU02x,1/upU02x,upU02mass,upU02green,upU02cost);
end

function upU02value=upU02heat(upU02t,upU02x)
if upU02t==0
    upU02value=0;
elseif upU02t<0.12
    upU02m=-4:4;
    upU02value=sum(exp(-(upU02x-4*upU02m).^2/(4*upU02t))-exp(-(upU02x-(4*upU02m+2)).^2/(4*upU02t)))/sqrt(4*pi*upU02t);
else
    upU02odd=1:2:79;
    upU02value=sum(exp(-(pi*upU02odd/2).^2*upU02t).*cos(pi*upU02odd*upU02x/2));
end
end

function upU02value=upU02survival(upU02t)
if upU02t==0
    upU02value=1;
elseif upU02t<0.12
    upU02k=0:4;
    upU02value=1-2*sum((-1).^upU02k.*erfc((2*upU02k+1)/(2*sqrt(upU02t))));
else
    upU02k=0:39;
    upU02odd=2*upU02k+1;
    upU02value=4/pi*sum((-1).^upU02k./upU02odd.*exp(-(pi*upU02odd/2).^2*upU02t));
end
end
