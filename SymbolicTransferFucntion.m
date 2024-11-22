syms m1 m2 m3 s x1 x2 x3 k1 k2 k3 c1 c2 c3 F1prime F2 
F1prime == (m1*s^2 + (k1 + k2) + (c1 + c2)*s )*x1 - (k2 + c2*s)*x2
F2 == (m2*s^2 + (k2 + k3) + (c2+c3)*s)*x2 - (k2 + c2*s)*x1 - (k3 + c3 * s) * x3
0 == (m3*s^2 + k3 + c3 * s) * x3 - ( k3 + c3*s)*x2

G11 = solve(x1/F1prime)
G11