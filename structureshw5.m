syms a b c y1 y2 y3 p q
E1 = y1*(b^3/3 + b^2*a) + y2*(b^2*c/2 + b*c*a - b*a^2/2) - y3*(b^2/2 + b*a) + p*b*a^2/2 -a*b*q == 0;
E2 = y2*(c^3/3 + c^2*b - (c - a)^3/3) + y3*(-c^2/2 - c*b - c*a + a^2/2) + y1*(c*b^2/2 + b*a*c - b*a^2/2) + p*(c*a^2/2 - a^3/3) -q*(c*a - a^2/2) == 0;
E3 = y2*(c^2/2 + b*c + c*a - a^2/2) - y3*(c + a + b) + y1*(b^2/2 + b*a) + p*a^2/2 - q*a == 0;

result = solve(E1,E2,E3);

% Y1, Y2, Y3 as a function of q as well
Y1_ans = vpa(result.y1);
Y2_ans = vpa(result.y2);
Y3_ans = vpa(result.y3);

% Y1, Y2, Y3 values when q = 0 is substituted 
Y1_final = subs(Y1_ans,q,0)
Y2_final = subs(Y2_ans,q,0)
Y3_final = subs(Y3_ans,q,0)

syms s3 s1 s2 E I
M1 = Y2_ans*s1 - Y3_ans;
M2 = Y1_ans*s2 + Y2_ans*c - Y3_ans;
M3 = p*s3 + Y1_ans*b + Y2_ans*(c - s3) - Y3_ans - q;
dM1 = diff(M1,q); dM2 = diff(M2,q); dM3 = diff(M3,q);
coeff1 = int(M1*dM1,s1,0,c);
coeff2 = int(M2*dM2,s2,0,b);
coeff3 = int(M3*dM3,s3,0,a);

% final theta value when q = 0 is substituted 
theta_final = simplify(1/E/I*subs(coeff1 + coeff2 + coeff3,q,0))

