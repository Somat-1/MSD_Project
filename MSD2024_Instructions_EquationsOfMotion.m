%% MSD 2024 Assignment Part 1: Deriving Transfer Functions using Equations of Motion
% Script prepared by: Aditya Natu, PhD Candidate, Mechatronics System Design

%% 
clear
clc

% Define bode options
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
opts.XLim= [1 1000];
opts.PhaseVisible = 'on';
opts.PhaseWrapping = 'on';
opts.PhaseMatching = 'on';
opts.PhaseMatchingFreq = 0;
opts.Grid = 'on';
opts.XLabel.String = 'Frequency';
opts.XLabel.FontSize = 16;
opts.XLabel.FontWeight = 'bold';
opts.YLabel.String = {'Magnitude','Phase'};
opts.YLabel.FontSize = 16;
opts.YLabel.FontWeight = 'bold';
opts.TickLabel.FontSize = 16;

%% Symbollically compute Equations of Motion (EOM)
% Create s
syms s x1 x2 x3 dx1 dx2 dx3 ddx1 ddx2 ddx3 m1 m2 m3 c1 c2 c3 k1 k2 k3 F1 F2 F3

% Define System Matrices (based on FBD)
% Rearrange in General Form: M*Xdd + C*Xd + K*X = F

% Mass and Stiffness matrix
M = [m1 0 0; 0 m2 0; 0 0 m3];
C = [c1+c2 -c2 0; -c2 c2+c3 -c3; 0 -c3 c3]; % Damping can be ignored
K = [k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3];
F = [F1 - F2; F2; 0];  % External forces (Only 2 Actuation Forces)

X = [x1; x2; x3];  % Displacement
Xd = [dx1; dx2; dx3];  % Velocity
Xdd = [ddx1; ddx2; ddx3];  % Acceleration

% General Form of Equations of Motion (EOM)
EOM = M*Xdd + C*Xd + K*X - F;
EOM = subs(EOM, [dx1 dx2 dx3 ddx1 ddx2 ddx3], [s*x1 s*x2 s*x3 s^2*x1 s^2*x2 s^2*x3]);

%% Defining Parameters for the transfer function x/F1 (wrt 1st Actuator Force)
c = 0.1e0;
EOMtf = subs(EOM, [m1 m2 m3 c1 c2 c3 k1 k2 k3], [2 0.2 0.05 c c c 1e4 3e4 4e4]);
EOMtf = subs(EOMtf, F2, 0);  % F2 is set to 0

% The system is a set of equations for x1, x2, and x3.
% Create the system of equations
eqs = EOMtf;

% Solve for x1, x2, and x3 in terms of F1
sol = solve(eqs, [x1, x2, x3]);

% Extract x in terms of F1
x1_sol = sol.x1;
x2_sol = sol.x2;
x3_sol = sol.x3;

% Now, compute the transfer function H(s) = x1/F1
H_x1F1 = simplify(x1_sol / F1);
H_x2F1 = simplify(x2_sol / F1);
H_x3F1 = simplify(x3_sol / F1);

% Display the transfer function
disp('Transfer Function x/F1:');
disp(H_x1F1);

% Extract the coefficients of the numerator and denominator:
[num, den] = numden(simplify(x1_sol / F1));
[num2, den2] = numden(simplify(x2_sol / F1));
[num3, den3] = numden(simplify(x3_sol / F1));

% Now, you can create the transfer function using the coefficients:
H_x1F1 = tf([sym2poly(num)], [sym2poly(den)]);
H_x2F1 = tf([sym2poly(num2)], [sym2poly(den2)]);
H_x3F1 = tf([sym2poly(num3)], [sym2poly(den3)]);

% Display the transfer function
disp('Transfer Function x1/F1:');
disp(H_x1F1);

figure(1);clf(1);
bode(H_x1F1,H_x2F1,H_x3F1,opts)
title('Bode plot EOM, x1/F1 ')
legend

%% Defining Parameters for the transfer function x/F2 (wrt 2nd Actuator Force)
c = 0.1e0;
EOMtf = subs(EOM, [m1 m2 m3 c1 c2 c3 k1 k2 k3], [2 0.2 0.05 c c c 1e4 3e4 4e4]);
EOMtf = subs(EOMtf, F1, 0);  % F1 is set to 0 

% Create the system of equations
eqs = EOMtf;

% Solve for x1, x2, and x3 in terms of F2
sol = solve(eqs, [x1, x2, x3]);

% Extract x1 in terms of F1 (assuming we want x1 only)
x1_sol = sol.x1;
x2_sol = sol.x2;
x3_sol = sol.x3;

% Now, compute the transfer function H(s) = x1/F1
H_x1F2 = simplify(x1_sol / F2);
H_x2F2 = simplify(x2_sol / F2);
H_x3F2 = simplify(x3_sol / F2);

% Extract the coefficients of the numerator and denominator:
[num, den] = numden(simplify(x1_sol / F2));
[num2, den2] = numden(simplify(x2_sol / F2));
[num3, den3] = numden(simplify(x3_sol / F2));

% Now, you can create the transfer function using the coefficients:
H_x1F2 = tf([sym2poly(num)], [sym2poly(den)]);
H_x2F2 = tf([sym2poly(num2)], [sym2poly(den2)]);
H_x3F2 = tf([sym2poly(num3)], [sym2poly(den3)]);

figure(2);clf(2);
bode(H_x1F2,H_x2F2,H_x3F2,opts)
title('Bode plot EOM x/F2')
legend

%% Create System Transfer Function Matrix

%[x] = [H]*[F] :- Relatting H to input and output vectors
H = [H_x1F1, H_x1F2;
     H_x2F1, H_x2F2;
     H_x3F1, H_x3F2];

% Input-Output Names
H.u{1} = 'F_1'; H.u{2} = 'F_2';
H.y{1} = 'x_1'; H.y{2} = 'x_2'; H.y{3} = 'x_3';

% Plot the System Matrix
figure(3);clf(3);
bode(H,opts)

% Customize figure appearance
set(gca, 'Color', 'w');
fig = gcf; 
fig.Color = 'w';
% Maximize the figure window
set(fig, 'WindowState', 'maximized');