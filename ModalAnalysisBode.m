clear; clc;

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

% Mass and stiffness matrices
M = [m1 0 0;
     0 m2 0;
     0 0 m3];
K = [(k1+k2) -k2 0;
    -k2 (k2+k3) -k3;
    0 -k3 k3];

% Solve for eigenvalues (T) and eigenvectors (phi)
[phi, T] = eig(-M^-1 * K);
phi1 = phi(:,1);
phi2 = phi(:,2);
phi3 = phi(:,3);

MM1 = phi1'*M*phi1;
MM2 = phi2'*M*phi2;
MM3 = phi3'*M*phi3;

MK1 = phi1'*K*phi1;
MK2 = phi2'*K*phi2;
MK3 = phi3'*K*phi3;

% Eigenmodes
f01 = sqrt(MM1\MK1)/(2*pi);
f02 = sqrt(MM2\MK2)/(2*pi);
f03 = sqrt(MM3\MK3)/(2*pi);

% Display results
disp('Eigenvector matrix (phi):');
disp(phi);
disp('Modal Masses:');
disp([MM1, MM2, MM3]);
disp('Modal Stiffness:');
disp([MK1, MK2, MK3]);
disp('Eigenmodes (Hz):');
disp([f01, f02, f03]);

% Define Laplace variable
s = tf('s');

% Modal Masses and Stiffnesses
MM = [MM1, MM2, MM3];
MK = [MK1, MK2, MK3];

% Modal Participation Factors
phiMatrix = phi;

% Pre-allocate transfer function matrix [l by k]
TransferFunctionMatrix = cell(3, 2);

% Generate transfer functions for l = [1,2,3], k = [1,2]
for l = 1:3
    for k = 1:2
        % Initialize transfer function
        totalTransferFunction = 0;
        
        % Summation over i = 1:3
        for i = 1:3
            % Extract modal properties for mode i
            phi_i_k = phiMatrix(k, i); % k-th element of i-th eigenvector
            phi_i_l = phiMatrix(l, i); % l-th element of i-th eigenvector
            M_i = MM(i);              % Modal mass
            K_i = MK(i);              % Modal stiffness
            
            % Transfer function term
            term = (phi_i_k * phi_i_l) / (M_i * s^2 + K_i);
            
            % Add to total transfer function
            totalTransferFunction = totalTransferFunction + term;
        end
        
        % Store in the transfer function matrix
        TransferFunctionMatrix{l, k} = totalTransferFunction;
    end
end

% Display the transfer function matrix
disp('Transfer Function Matrix [l by k]:');
for l = 1:3
    for k = 1:2
        fprintf('T(%d,%d):\n', l, k);
        disp(TransferFunctionMatrix{l, k});
    end
end

% Plot Bode diagrams for each transfer function
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz';
w = logspace(-1, 3, 1001) * 2 * pi; % Frequency range

figure;
plotIdx = 1; % Subplot index

for l = 1:3
    for k = 1:2
        % Extract transfer function
        TF = TransferFunctionMatrix{l, k};
        
        % Create individual Bode plots
        subplot(3, 2, plotIdx);
        bode(TF, w, opts);
        grid on;
        title(sprintf('Bode Plot for T(%d,%d)', l, k));
        plotIdx = plotIdx + 1; % Increment subplot index
    end
end
