Fs = 8192;
T  = 10;
N  = Fs * T;             
t  = (0:N-1)' / Fs;     

Vload = out.simout.Data(:,1);
Vshift = circshift(Vload,-7);

LoD = [ ...
   -0.07576571478927333, ...
   -0.02963552764599851, ...
    0.49761866763201545, ...
    0.80373875180591610, ...
    0.29785779560527736, ...
   -0.09921954357684722, ...
   -0.01260396726203783, ...
    0.03222310060404270 ];

HiD = [ ...
   -0.03222310060404270, ...
   -0.01260396726203783, ...
    0.09921954357684722, ...
    0.29785779560527736, ...
   -0.80373875180591610, ...
    0.49761866763201545, ...
    0.02963552764599851, ...
   -0.07576571478927333 ];


% ============================================================
%       ASSIGNMENT 5.3: Square Wave 6-LEVEL DWT (sym4)
% ============================================================

max_level = 6;

% Cell arrays to hold approximation and detail vectors across all levels
% A_coeff{1} = Vshift (Level 0), A_coeff{7} = A6
% D_coeff{1} = D1, ..., D_coeff{6} = D6
A_coeff = cell(max_level + 1, 1);
D_coeff = cell(max_level, 1);
A_coeff{1} = Vload;

% ------------------ DECONSTRUCTION LOOP ------------------
for j = 1:max_level
    A_input_shifted = circshift(A_coeff{j}, -7);
    N_prev = length(A_input_shifted);
    N_curr = N_prev / 2;
    
    A_temp = zeros(N_curr, 1);
    D_temp = zeros(N_curr, 1);
    
    for n = 0:(N_curr - 1)
        sumA = 0;
        sumD = 0;
        for k = 0:7
            index = mod(2*n + k, N_prev) + 1; % Circular boundary wrapping
            sumA = sumA + LoD(k+1) * A_input_shifted(index);
            sumD = sumD + HiD(k+1) * A_input_shifted(index);
        end
        A_temp(n+1) = sumA;
        D_temp(n+1) = sumD;
    end
    
    A_coeff{j+1} = A_temp;
    D_coeff{j}   = D_temp;
end

% Extract final Level 6 coefficients (Length = 81920 / 2^6 = 1280 samples)
% Store target coefficients and labels in cell arrays
coeffs = {A_coeff{7}, D_coeff{6}, D_coeff{5}, D_coeff{4}, D_coeff{3}, D_coeff{2}, D_coeff{1}};
labels = {'A6', 'D6', 'D5', 'D4', 'D3', 'D2', 'D1'};

fprintf('============================================\n');
fprintf('     6-LEVEL DECOMPOSITION CONFIRMATION     \n');
fprintf('============================================\n');
fprintf('Length of Vload_shifted : %d samples\n', length(Vshift));
fprintf('Length of D6            : %d coefficients\n', length(D_coeff{1}));
fprintf('Length of D5            : %d coefficients\n', length(D_coeff{2}));
fprintf('Length of D4            : %d coefficients\n', length(D_coeff{3}));
fprintf('Length of D3            : %d coefficients\n', length(D_coeff{4}));
fprintf('Length of D2            : %d coefficients\n', length(D_coeff{5}));
fprintf('Length of D1            : %d coefficients\n', length(D_coeff{6}));
fprintf('Length of A6            : %d coefficients\n', length(A_coeff{7}));


% Export A6 and D6 to CSV
A6_D6_Table = table(A_coeff{7}, D_coeff{6}, 'VariableNames', {'A6', 'D6'});
writetable(A6_D6_Table, 'A6_D6_Coefficients.csv');

% ------------------ RECONSTRUCTION LOOP ------------------
A_rec = A_coeff{7};
for j = max_level:-1:1
    N_curr = length(A_rec);
    N_prev = N_curr * 2;
    A_prev_rec = zeros(N_prev, 1);
    
    for n = 0:(N_curr - 1)
        for k = 0:7
            index = mod(2*n + k + 7, N_prev) + 1;
            A_prev_rec(index) = A_prev_rec(index) + ...
                A_rec(n+1) * LoD(k+1) + D_coeff{j}(n+1) * HiD(k+1);
        end
    end
    A_rec = A_prev_rec;
end

Vload_rec   = A_rec;
recon_error = max(abs(Vload - Vload_rec));

fprintf('\n============================================\n');
fprintf('Maximum Reconstruction Error (Level 6): %e\n', recon_error);
fprintf('============================================\n');

% =========================================================
%      ENERGY AND PERCENTAGE DISTRIBUTION CALCULATIONS
% =========================================================

% 1. Calculate Total Signal Energy in time-domain
total_energy = sum(Vload.^2);

energies     = zeros(1, 7);
pct_energies = zeros(1, 7);

% 2. Calculate and display formatted results
fprintf('\n=========================================================\n');
fprintf('     COEFFICIENT ENERGY AND PERCENTAGE DISTRIBUTION      \n');
fprintf('=========================================================\n');
fprintf('%-6s | %-12s | %-18s | %-12s\n', 'Coeff', 'Length', 'Energy (V^2)', '% Energy');
fprintf('---------------------------------------------------------\n');

for i = 1:length(coeffs)
    energies(i)     = sum(coeffs{i}.^2);                   % Sum of squared coefficients
    pct_energies(i) = (energies(i) / total_energy) * 100;  % Percentage of total energy
    
    fprintf('%-6s | %-12d | %-18.4e | %-11.4f%%\n', ...
        labels{i}, length(coeffs{i}), energies(i), pct_energies(i));
end

fprintf('---------------------------------------------------------\n');
fprintf('=========================================================\n');

% ------------------ PLOTTING RESULTS ------------------
figure('Name', '6-Level Wavelet Decomposition', 'NumberTitle', 'off');

% Plot Vshift
subplot(2,1,1);
plot(t, Vshift, 'b', 'LineWidth', 1);
title('Load Voltage (Measured Signal)');
xlabel('Time (s)'); ylabel('Voltage (V)');
grid on;

% Downsampled time vector for Level 6 (1280 samples)

% Plot A6
subplot(2,1,2);
plot(t, Vload_rec, 'r', 'LineWidth', 1);
title('Level 6 Approximation Reconstructed Signal');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

