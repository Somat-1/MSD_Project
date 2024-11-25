clear; clc;

%define symbols 
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

%%%to neglect damping 
%c1 = 0;   % Ns/m
%c2 = 0;   % Ns/m
%c3 = 0;   % Ns/m

% define and solve equations of motion 
eq1 = 0 == (m1*s^2 + (k1 + k2) + (c1 + c2)*s )*x1 - (k2 + c2*s)*x2 - F1 + F2; %EOM block 1
eq2 = F2 == (m2*s^2 + (k2 + k3) + (c2+c3)*s)*x2 - (k2 + c2*s)*x1 - (k3 + c3 * s) * x3;  %EOM block 2
eq3 = 0 == (m3*s^2 + k3 + c3 * s) * x3 - ( k3 + c3*s)*x2;  %EOM block 3

S = solve(eq1,eq2,eq3); % Solves EOM for x1, x2, x3 in terms of F1, F2 

%%% call x1, x2, x3 solutions  
expr_x1 = simplify(S.x1);
expr_x2 = simplify(S.x2);
expr_x3 = simplify(S.x3);

%%% Go from x to transfer function 
G11 = x_to_G(expr_x1, F1, F2)
G21 = x_to_G(expr_x2, F1, F2)
G31 = x_to_G(expr_x3, F1, F2)
G12 = x_to_G(expr_x1, F2, F1)
G22 = x_to_G(expr_x2, F2, F1)
G32 = x_to_G(expr_x3, F2, F1)

%bode plot of x1/F1
subplot(3,2,1);
bode(G11)
grid on
xlabel('Frequency (rad/sec)');
ylabel('Magnitude (dB)');
title('Bode Diagram G11');
hold on;

%bode plot of x2/F1
subplot(3,2,3);
bode(G21)
grid on
xlabel('Frequency (rad/sec)');
ylabel('Magnitude (dB)');
title('Bode Diagram G21');
hold on;

%bode plot of x3/F1
subplot(3,2,5);
bode(G31)
grid on
xlabel('Frequency (rad/sec)');
ylabel('Magnitude (dB)');
title('Bode Diagram G31');
hold on;

%bode plot of x1/F2
subplot(3,2,2);
bode(G12)
grid on
xlabel('Frequency (rad/sec)');
ylabel('Magnitude (dB)');
title('Bode Diagram G12');
hold on;

%bode plot of x2/F2
subplot(3,2,4);
bode(G22)
grid on
xlabel('Frequency (rad/sec)');
ylabel('Magnitude (dB)');
title('Bode Diagram G22');
hold on;

%bode plot of x3/F2
subplot(3,2,6);
bode(G32)
grid on
xlabel('Frequency (rad/sec)');
ylabel('Magnitude (dB)');
title('Bode Diagram G32');
hold on;













