% Primary serialized only, independent n4r11 prefix, no clear/path.
n4r11P=[1+sqrt(2)/2,(2+2^(3/4)*3^(1/4))/4,(1+2^(1/4))/2,(1+2^(1/4))/2];
n4r11Q=n4r11P.*[2 2 1 1]; n4r11h=[.01 .02 .01 0]; n4r11d=n4r11P.*n4r11h;
n4r11base=[0 0 0 0;0 .5 0 0;.5 0 0 0;.5 .5 0 0]; n4r11sign=[1;-1;-1;1];
n4r11u=n4r11h./[2 2 1 1]; n4r11phase=n4r11base+n4r11sign*n4r11u;
n4r11physical=n4r11base.*n4r11Q+n4r11sign*n4r11d;
assert(max(max(abs(n4r11phase.*n4r11Q-n4r11physical)))<1e-12);
n4r11counts=[24 16 36 24]; n4r11eps=.05; n4r11weights=ones(4,1)+n4r11eps*n4r11sign;
n4r11mass=n4r11counts*n4r11weights; n4r11r=n4r11mass^(1/4)/2;
assert(abs(n4r11mass-99.8)<1e-12 && sum(n4r11weights)==4);
fprintf('n4 round011 counts=[%d %d %d %d],delta=%d,mass=%.16g,r=%.16g,cellmass=%.16g\n',n4r11counts,n4r11counts*n4r11sign,n4r11mass,n4r11r,sum(n4r11weights));
