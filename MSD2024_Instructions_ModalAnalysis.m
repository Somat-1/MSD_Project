%% MSD 2024 Assignment Part 1: Modal Analysis
% Script prepared by: Aditya Natu, PhD Candidate, Mechatronics System Design

%% Modal Analysis 
%clear all

% System Parameters
m1 = 2;    % [kg] mass of the base
m2 = 0.2;  % [kg] mass of the manipulator
m3 = 0.05; % [kg] mass of the parastic part

k1 = 1e4;  % [N/m] stiffness coefficient between the ground and the base
k2 = 3e4;  % [N/m] stiffness coefficient between the base and the manipulator
k3 = 4e4;  % [N/m] stiffness coefficient between the manipulator and parasitic part

c1 = 0.1e0; % [Ns/m] damping coefficient between the ground and the base
c2 = 0.1e0; % [Ns/m] damping coefficient between the base and the manipulator
c3 = 0.1e0; % [Ns/m] damping coefficient between the manipulator and parasitic part

% Mass, Damping and Stiffness matrix
M = [m1 0 0; 0 m2 0; 0 0 m3];
C = [c1+c2 -c2 0;-c2 c2+c3 -c3;0 -c3 c3]; % Damping can be ignored
K = [k1+k2 -k2 0;-k2 k2+k3 -k3;0 -k3 k3];

% Eigenvalues and Eigenvector computation
[phi, D]=eig(M^-1*K);

phi1=phi(:,1); %eigenvector of mode 1
phi2=phi(:,2); %eigenvector of mode 2
phi3=phi(:,3); %eigenvector of mode 3

MM1=phi1'*M*phi1; % Modal Mass 1
MM2=phi2'*M*phi2; % Modal Mass 2
MM3=phi3'*M*phi3; % Modal Mass 3

MK1=phi1'*K*phi1; % Modal Stiffness 1
MK2=phi2'*K*phi2; % Modal Stiffness 2
MK3=phi3'*K*phi3; % Modal Stiffness 3

MC1=phi1'*C*phi1; % Modal Damping 1
MC2=phi2'*C*phi2; % Modal Damping 2
MC3=phi3'*C*phi3; % Modal Damping 3

% Eigenfrequencies
f01=sqrt(MM1\MK1)/(2*pi); % Eigenfrequency 1
f02=sqrt(MM2\MK2)/(2*pi); % Eigenfrequency 2
f03=sqrt(MM3\MK3)/(2*pi); % Eigenfrequency 3


% Define bode options
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
opts.XLim= [1 1000];
opts.PhaseVisible = 'on';
opts.PhaseWrapping = 'on';
opts.Grid = 'on';
opts.XLabel.String = 'Frequency';
opts.XLabel.FontSize = 16;
opts.XLabel.FontWeight = 'bold';
opts.YLabel.String = {'Magnitude','Phase'};
opts.YLabel.FontSize = 16;
opts.YLabel.FontWeight = 'bold';
opts.TickLabel.FontSize = 16;

s=tf('s'); w=logspace(-1,4,20000)*2*pi;

%% Transfer Functions from Modal Analysis
% General Form
% xl/Fk (for i^th mode) = phi_i(l)*phi_i(k)/(MM_i*s^2 + MK_i)

% Transfer Function of x1/F1 (l=1, k=1)
xF11=phi1(1)*phi1(1)/(MM1*s^2+MK1)+phi2(1)*phi2(1)/(MM2*s^2+MK2)+phi3(1)*phi3(1)/(MM3*s^2+MK3) % Damping Excluded
xF11d=phi1(1)*phi1(1)/(MM1*s^2+MC1*s+MK1)+phi2(1)*phi2(1)/(MM2*s^2+MC2*s+MK2)+phi3(1)*phi3(1)/(MM3*s^2+MC3*s+MK3); % Damping Included
%close all
figure; clf(1); bode(xF11,w,opts); grid on
hold on 

%% Plots of individual modes and TF 
xF11m1=phi1(1)*phi1(1)/(MM1*s^2+MK1);
bode(xF11m1,w,opts); grid on; hold on;
xF11m2=phi2(1)*phi2(1)/(MM2*s^2+MK2);
bode(xF11m2,w,opts); grid on; hold on;
xF11m3=phi3(1)*phi3(1)/(MM3*s^2+MK3);
%bode(xF11m3,w,opts); grid on;
title('x1F1 and the 3 modes it comprises of')
%hold on; bode(xF11d,w,opts); grid on; % uncomment to add undamped one to
%plot
%title('Bode F1/x1 damped and undamped')

legend

%% Other DOF wrt F1
% Transfer Function of x2/F1 (l=2, k=1)
xF21=phi1(1)*phi1(2)/(MM1*s^2+MK1)+phi2(1)*phi2(2)/(MM2*s^2+MK2)+phi3(1)*phi3(2)/(MM3*s^2+MK3)
figure(1); clf(1); bode(xF21,w,opts); grid on; hold on;

xF21m1=phi1(1)*phi1(2)/(MM1*s^2+MK1);
bode(xF21m1,w,opts); grid on; hold on;
xF21m2=phi2(1)*phi2(2)/(MM2*s^2+MK2);
bode(xF21m2,w,opts); grid on; hold on;
xF21m3=phi3(1)*phi3(2)/(MM3*s^2+MK3);
bode(xF21m3,w,opts); grid on;
legend
title('x2F1 and the 3 modes it comprises of')
% Transfer Function of x3/F1 (l=3, k=1)
xF31=phi1(1)*phi1(3)/(MM1*s^2+MK1)+phi2(1)*phi2(3)/(MM2*s^2+MK2)+phi3(1)*phi3(3)/(MM3*s^2+MK3);
bode(xF31,w,opts); grid on
legend

%% Transfer Functions wrt F2
% Transfer Function of x1/F2 (l=1, k=2)
l=1;k=2;
xF12=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);

xF12f = xF12 - xF11; % Net Force on m1 is difference of two actuator forces (F1-F2)
figure(2); clf(2); 
bode(xF12f,w,opts); grid on; hold on;

% Transfer Function of x2/F2 (l=2, k=2)
l=2;k=2;
xF22=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);
xF22f = xF22 - xF21;
bode(xF22f,w,opts); grid on

% Transfer Function of x3/F2 (l=3, k=2)
l=3;k=2;
xF32=phi1(l)*phi1(k)/(MM1*s^2+MK1)+phi2(l)*phi2(k)/(MM2*s^2+MK2)+phi3(l)*phi3(k)/(MM3*s^2+MK3);
xF32f = xF32 - xF31;
bode(xF32f,w,opts); grid on


%% %% Create System Transfer Function Matrix

%[x] = [MA1]*[Fnet] :- Relatting MA1 to input (Actuation Forces on Body) and output vectors
MA1 = [xF11, xF12;
       xF21, xF22;
       xF31, xF32];

%[x] = [MA2]*[F] :- Relatting  MA2to input (Actuator Forces) and output vectors
MA2 = [xF11, xF12f;
       xF21, xF22f;
       xF31, xF32f];

opts.PhaseVisible = 'off';

figure(3);clf(3);
bode(MA1,'k',MA2,'r--',opts)
%% Compare Matrices from EOM and Modal Analysis

figure(4);clf(4);
bode(H,'k',MA2,'r:',opts)
legend Location southwest
title('eom vs modal')
% Customize figure appearance
set(gca, 'Color', 'w');
fig = gcf; 
fig.Color = 'w';
% Maximize the figure window
set(fig, 'WindowState', 'maximized');