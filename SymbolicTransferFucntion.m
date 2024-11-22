clear; clc;

%define symbols 
syms m1 m2 m3 s x1 x2 x3 k1 k2 k3 c1 c2 c3 F1prime F2 

% define equations of motion 
eq1 = F1prime == (m1*s^2 + (k1 + k2) + (c1 + c2)*s )*x1 - (k2 + c2*s)*x2; %EOM block 1
eq2 = F2 == (m2*s^2 + (k2 + k3) + (c2+c3)*s)*x2 - (k2 + c2*s)*x1 - (k3 + c3 * s) * x3;  %EOM block 2
eq3 = 0 == (m3*s^2 + k3 + c3 * s) * x3 - ( k3 + c3*s)*x2;  %EOM block 3

S = solve(eq1,eq2,eq3);

S.x1
S.F1prime

%H11 = simplify(S.x1/F1prime)












