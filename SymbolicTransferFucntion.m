clear; clc;

%define symbols 
%syms m1 m2 m3 s x1 x2 x3 k1 k2 k3 c1 c2 c3 F1 F2 
syms F1 F2 s x1 x2 x3 

% Masses
m1 = 2;     % kg
m2 = 0.2;   % kg
m3 = 0.05;  % kg

% Spring constants
k1 = 10000; % N/m
k2 = 30000; % N/m
k3 = 40000; % N/m

% Damping coefficients
c1 = 0.1;   % Ns/m
c2 = 0.1;   % Ns/m
c3 = 0.1;   % Ns/m

%to neglect damping 
%c1 = 0;   % Ns/m
%c2 = 0;   % Ns/m
%c3 = 0;   % Ns/m

% define and solve equations of motion 
eq1 = 0 == (m1*s^2 + (k1 + k2) + (c1 + c2)*s )*x1 - (k2 + c2*s)*x2 - F1 + F2; %EOM block 1
eq2 = F2 == (m2*s^2 + (k2 + k3) + (c2+c3)*s)*x2 - (k2 + c2*s)*x1 - (k3 + c3 * s) * x3;  %EOM block 2
eq3 = 0 == (m3*s^2 + k3 + c3 * s) * x3 - ( k3 + c3*s)*x2;  %EOM block 3

S = solve(eq1,eq2,eq3); % Solves EOM for x1, x2, x3 in terms of F1, F2 


%%% Split the solutions for x1, x2, x3, in a solution ito F1 + F2 aka G11 and G12  
expr = simplify(S.x1);

% First expand the expression
expr_expanded = expand(expr);

% Collect terms with F1
F1_terms = collect(expr_expanded, F1);

% Collect terms with F2
F2_terms = expr_expanded - F1_terms;

% Or alternatively, you can use these commands:
F1_part = subs(expr_expanded, F2, 0);  % This will give you only F1 terms
F2_part = subs(expr_expanded, F1, 0); % This will give you only F2 terms

simplify(F1_part);
simplify(F2_part);

F1_part = F1_part/F1
F2_part = F2_part/F2

%%% split G11 in numenator and denom and make arrays with coefficients 
[num, den] = numden(F1_part);

% Get coefficients of numerator
[num_coeffs, num_terms] = coeffs(num, s);

% Get coefficients of denominator
[den_coeffs, den_terms] = coeffs(den, s);

% Convert symbolic coefficients to numeric double precision
num_coeffs_numeric = double(num_coeffs);
den_coeffs_numeric = double(den_coeffs);

% Ensure they are row vectors
if ~isrow(num_coeffs_numeric)
    num_coeffs_numeric = num_coeffs_numeric';
end
if ~isrow(den_coeffs_numeric)
    den_coeffs_numeric = den_coeffs_numeric';
end

% Create transfer function
F1_tf = tf(num_coeffs_numeric, den_coeffs_numeric)

%bode plot of x1/F1
bode(F1_tf)
grid










