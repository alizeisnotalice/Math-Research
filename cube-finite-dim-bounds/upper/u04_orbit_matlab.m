% U04: independent actual Delta Green orbit check, no workspace reset.
upU04R=2.1;
upU04theta=pi/8;
upU04c=cos(upU04theta);
upU04s=sin(upU04theta);
upU04points=[1,1,1;1,0,0;1,1,0;1,.7,.3];
upU04mass=integral(@(upU04t) arrayfun(@(upU04z) upU04survival(upU04z,upU04R)^3,upU04t),0,Inf,'AbsTol',1e-12,'RelTol',1e-10);
for upU04q=1:size(upU04points,1)
    upU04point=upU04points(upU04q,:);
    upU04value=0;
    for upU04i=1:3
        for upU04j=1:3
            if upU04i==upU04j
                continue
            end
            upU04k=6-upU04i-upU04j;
            for upU04sign=[-1,1]
                upU04y=[upU04c*upU04point(upU04i)+upU04s*upU04sign*upU04point(upU04j), ...
                    -upU04s*upU04point(upU04i)+upU04c*upU04sign*upU04point(upU04j),upU04point(upU04k)];
                upU04value=upU04value+integral(@(upU04t) arrayfun(@(upU04z) prod(arrayfun(@(upU04coord) upU04heat(upU04z,upU04coord,upU04R),upU04y)),upU04t),0,Inf,'AbsTol',1e-12,'RelTol',1e-10)/12;
            end
        end
    end
    fprintf('U04 point=(%.3g,%.3g,%.3g) orbitGreen=%.17g mass=%.17g pointCost=%.17g\n',upU04point,upU04value,upU04mass,upU04mass/(8*upU04value));
end
upU04B3=(1/sqrt(3)-1/upU04R)/(4*pi);
upU04B4=(1/4-1/upU04R^2)/(4*pi^2);
fprintf('U04 global cube lower B3=%.17g B4=%.17g coarse mass-bound costs=%.17g %.17g\n',upU04B3,upU04B4,upU04R^2/(16*upU04B3),upU04R^2/(32*upU04B4));

function upU04v=upU04heat(upU04t,upU04x,upU04R)
if upU04t==0 || abs(upU04x)>=upU04R
    upU04v=0;
elseif upU04t/upU04R^2<.12
    upU04m=-4:4;
    upU04v=sum(exp(-(upU04x-4*upU04m*upU04R).^2/(4*upU04t))-exp(-(upU04x-(4*upU04m+2)*upU04R).^2/(4*upU04t)))/sqrt(4*pi*upU04t);
else
    upU04odd=1:2:79;
    upU04v=sum(exp(-(pi*upU04odd/(2*upU04R)).^2*upU04t).*cos(pi*upU04odd*upU04x/(2*upU04R)))/upU04R;
end
end

function upU04v=upU04survival(upU04t,upU04R)
if upU04t==0
    upU04v=1;
elseif upU04t/upU04R^2<.12
    upU04k=0:4;
    upU04v=1-2*sum((-1).^upU04k.*erfc((2*upU04k+1)*upU04R/(2*sqrt(upU04t))));
else
    upU04k=0:39;
    upU04odd=2*upU04k+1;
    upU04v=4/pi*sum((-1).^upU04k./upU04odd.*exp(-(pi*upU04odd/(2*upU04R)).^2*upU04t));
end
end
