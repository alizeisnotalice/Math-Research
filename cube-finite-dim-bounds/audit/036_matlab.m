% 036 finite independent algebra crosscheck; no path/clear or Green nodes.
aud036g=sym(2621323989)/sym(500000000000);
aud036I=sym(23129127)/sym(40000000000);
aud036U=sym(1979)/sym(250);
aud036c2=27*aud036I^2/aud036g^3;
aud036gap=aud036U^2-aud036c2;
assert(isAlways(aud036gap>0));
aud036weight=0;
for aud036i=0:127
 for aud036j=0:aud036i
  aud036weight=aud036weight+1+(aud036i~=aud036j);
 end
end
assert(aud036weight==16384);
aud036l=sym('aud036l','positive');aud036u=sym('aud036u','positive');
assert(isAlways(simplify((aud036u^2+aud036u*aud036l+aud036l^2)*(aud036u-aud036l)-(aud036u^3-aud036l^3))==0));
aud036atan5=sym(0);aud036atan239=sym(0);
for aud036k=0:12
 aud036atan5=aud036atan5+(-1)^aud036k*sym(1)/sym(5)^(2*aud036k+1)/(2*aud036k+1);
end
for aud036k=0:3
 aud036atan239=aud036atan239+(-1)^aud036k*sym(1)/sym(239)^(2*aud036k+1)/(2*aud036k+1);
end
assert(isAlways(16*aud036atan5-4*aud036atan239<sym(355)/113));
fprintf('036 PASS area1 chord identity Machin upper; cost2=%s gap=%s\n',char(aud036c2),char(aud036gap));
disp(vpa(sqrt(aud036c2),50));
