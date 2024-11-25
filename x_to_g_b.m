function [G_tf] = x_to_g_b(expr,F_used, F_notused1, F_notused2)
%%% This one does the same as the x_to_G functions but for the system
%%% withthe third force
syms s

% First expand the expression
expr_expanded = expand(expr);

% Splits solution for x into part depended on only 1 F
G_part = subs(expr_expanded, [F_notused1,F_notused2], [0,0]);  % This will give you only terms of the other F

simplify(G_part);

% Take F out of expression 
G_part = G_part/F_used;

%%% split G11 in numenator and denom and make arrays with coefficients 
[num, den] = numden(G_part);

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
G_tf = tf(num_coeffs_numeric, den_coeffs_numeric);

end