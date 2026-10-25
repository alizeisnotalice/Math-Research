u7audit_t = sym('u7audit_t','real');
u7audit_b = (2-u7audit_t)^3/6;
u7audit_out = sym(zeros(1,2));
for u7audit_n=3:4
 u7audit_out(u7audit_n-2)=100*(subs(diff(u7audit_b,u7audit_t,2),u7audit_t,sym(3)/2)*(sym(2)/3)^(u7audit_n-1)-2*(u7audit_n-1)*subs(u7audit_b,u7audit_t,sym(3)/2)*(sym(2)/3)^(u7audit_n-2));
end
disp(u7audit_out);
assert(isequal(u7audit_out,[sym(50)/3,sym(250)/27]));
disp('ROUND007_AUDIT_SYMBOLIC_PASS');
