Fs = 8192;
T  = 10;
N  = Fs * T;             
t  = (0:N-1)' / Fs;     
Vload  = out.simout.Data(:,1);
Vshift = circshift(Vload,-7);
Iout   = out.simout.Data(:,2);
Ishift = circshift(Iout, -7);

LoD = [ ...
   -0.07576571478927333, -0.02963552764599851,  0.49761866763201545,  0.80373875180591610, ...
    0.29785779560527736, -0.09921954357684722, -0.01260396726203783,  0.03222310060404270 ];
HiD = [ ...
   -0.03222310060404270, -0.01260396726203783,  0.09921954357684722,  0.29785779560527736, ...
   -0.80373875180591610,  0.49761866763201545,  0.02963552764599851, -0.07576571478927333 ];

% ============================================================
%       DYNAMIC 200ms GATING SETUP (50 WINDOWS OVER 10s)
% ============================================================
Twin = 0.2;                         % 200 ms interval
num_windows = round(T / Twin);      % 50 windows total

q = cell(num_windows, 1);
for k = 1:num_windows
    g = zeros(N, 1);
    idx_start = round((k-1) * Fs * Twin) + 1;
    idx_end   = round(k * Fs * Twin);
    g(idx_start:idx_end) = 1;
    q{k} = g .* exp(-1i * 2 * pi * 50 * t);
end

% ============================================================
%       ASSIGNMENT 6: 6-LEVEL DWT (sym4) V & I
% ============================================================
max_level = 6;

% Cell arrays for Voltage
A_coeff_V = cell(max_level + 1, 1);
D_coeff_V = cell(max_level, 1);
A_coeff_V{1} = Vload;

% Cell arrays for Current
A_coeff_I = cell(max_level + 1, 1);
D_coeff_I = cell(max_level, 1);
A_coeff_I{1} = Iout;

% Preallocate cell arrays to hold the final H outputs
H_A6 = cell(num_windows, 1); 
H_D6 = cell(num_windows, 1);

% ------------------ DECONSTRUCTION LOOPS ------------------
for j = 1:max_level
    % Shift inputs
    A_input_shifted_V = circshift(A_coeff_V{j}, -7);
    A_input_shifted_I = circshift(A_coeff_I{j}, -7);
    
    N_prev = length(A_input_shifted_V);
    N_curr = N_prev / 2;
    
    A_temp_V = zeros(N_curr, 1);
    D_temp_V = zeros(N_curr, 1);
    A_temp_I = zeros(N_curr, 1);
    D_temp_I = zeros(N_curr, 1);
    
    for n = 0:(N_curr - 1)
        sumA_V = 0; sumD_V = 0;
        sumA_I = 0; sumD_I = 0;
        for k = 0:7
            index = mod(2*n + k, N_prev) + 1; 
            
            % Voltage Convolution
            sumA_V = sumA_V + LoD(k+1) * A_input_shifted_V(index);
            sumD_V = sumD_V + HiD(k+1) * A_input_shifted_V(index);
            
            % Current Convolution
            sumA_I = sumA_I + LoD(k+1) * A_input_shifted_I(index);
            sumD_I = sumD_I + HiD(k+1) * A_input_shifted_I(index);
        end
        A_temp_V(n+1) = sumA_V;
        D_temp_V(n+1) = sumD_V;
        A_temp_I(n+1) = sumA_I;
        D_temp_I(n+1) = sumD_I;
    end
    
    A_coeff_V{j+1} = A_temp_V;
    D_coeff_V{j}   = D_temp_V;
    A_coeff_I{j+1} = A_temp_I;
    D_coeff_I{j}   = D_temp_I;
end

% Deconstruction of "q"s to form "H"s for all 50 windows
for k = 1:num_windows
    A_coeff = cell(max_level + 1, 1);
    D_coeff = cell(max_level, 1);
    A_coeff{1} = q{k}; % Level 0 input
    
    for j = 1:max_level
        A_input_shifted = circshift(A_coeff{j}, -7);
        N_prev = length(A_input_shifted);
        N_curr = N_prev / 2;
        
        A_temp = zeros(N_curr, 1, 'like', 1i); % Complex vector allocation
        D_temp = zeros(N_curr, 1, 'like', 1i);
        
        for n = 0:(N_curr - 1)
            sumA = 0; sumD = 0;
            for idx_k = 0:7
                index = mod(2*n + idx_k, N_prev) + 1;
                sumA = sumA + LoD(idx_k+1) * A_input_shifted(index);
                sumD = sumD + HiD(idx_k+1) * A_input_shifted(index);
            end
            A_temp(n+1) = sumA;
            D_temp(n+1) = sumD;
        end
        
        A_coeff{j+1} = A_temp;
        D_coeff{j}   = D_temp;
    end
    
    % Save Level 6 Approximation (A6) and Detail (D6)
    H_A6{k} = A_coeff{7}; 
    H_D6{k} = D_coeff{6}; 
end

% Export A6 and D6 for Voltage and Current
A6_D6_Table = table(A_coeff_V{7}, D_coeff_V{6}, A_coeff_I{7}, D_coeff_I{6}, ...
    'VariableNames', {'Voltage_A6', 'Voltage_D6', 'Current_A6', 'Current_D6'});
writetable(A6_D6_Table, 'A6_D6_Coefficients.csv');

A5_Table = table(A_coeff_V{6}, A_coeff_I{6}, 'VariableNames', {'Voltage_A5', 'Current_A5'});
writetable(A5_Table, 'A5_Coefficients.xlsx');

% ------------------ RECONSTRUCTION LOOP ------------------
% Reconstruct Voltage
A_rec_V = A_coeff_V{7};
for j = max_level:-1:1
    N_curr = length(A_rec_V);
    N_prev = N_curr * 2;
    A_prev_rec_V = zeros(N_prev, 1);
    
    for n = 0:(N_curr - 1)
        for k = 0:7
            index = mod(2*n + k + 7, N_prev) + 1;
            A_prev_rec_V(index) = A_prev_rec_V(index) + ...
                A_rec_V(n+1) * LoD(k+1) + D_coeff_V{j}(n+1) * HiD(k+1);
        end
    end
    A_rec_V = A_prev_rec_V;
end
Vload_rec = A_rec_V;

% Reconstruct Current
A_rec_I = A_coeff_I{7};
for j = max_level:-1:1
    N_curr = length(A_rec_I);
    N_prev = N_curr * 2;
    A_prev_rec_I = zeros(N_prev, 1);
    
    for n = 0:(N_curr - 1)
        for k = 0:7
            index = mod(2*n + k + 7, N_prev) + 1;
            A_prev_rec_I(index) = A_prev_rec_I(index) + ...
                A_rec_I(n+1) * LoD(k+1) + D_coeff_I{j}(n+1) * HiD(k+1);
        end
    end
    A_rec_I = A_prev_rec_I;
end
Iout_rec = A_rec_I;

% Calculate Errors
recon_error_V = max(abs(Vload - Vload_rec));
recon_error_I = max(abs(Iout - Iout_rec));
fprintf('\n============================================\n');
fprintf('Maximum V Reconstruction Error (Level 6): %e\n', recon_error_V);
fprintf('Maximum I Reconstruction Error (Level 6): %e\n', recon_error_I);
fprintf('============================================\n');

% ============================================================
%      STAGE-BY-STAGE PHASOR AND IMPEDANCE CALCULATION
% ============================================================
CV = A_coeff_V{7}; % Voltage A6 Coefficients
CI = A_coeff_I{7}; % Current A6 Coefficients

% Preallocate arrays for 50 stages
V = zeros(num_windows, 1, 'like', 1i);
I = zeros(num_windows, 1, 'like', 1i);
Z = zeros(num_windows, 1, 'like', 1i);
R = zeros(num_windows, 1);
X = zeros(num_windows, 1);

% Compute 50 Hz Phasors and Impedance for each 200 ms window
for k = 1:num_windows
    V(k) = sum(CV .* H_A6{k});
    I(k) = sum(CI .* H_A6{k});
    Z(k) = V(k) / I(k);
    R(k) = real(Z(k));
    X(k) = imag(Z(k));
end

% =========================================================
%      ENERGY AND PERCENTAGE DISTRIBUTION CALCULATIONS
% =========================================================
coeffs_V = {A_coeff_V{7}, D_coeff_V{6}, D_coeff_V{5}, D_coeff_V{4}, D_coeff_V{3}, D_coeff_V{2}, D_coeff_V{1}};
coeffs_I = {A_coeff_I{7}, D_coeff_I{6}, D_coeff_I{5}, D_coeff_I{4}, D_coeff_I{3}, D_coeff_I{2}, D_coeff_I{1}};
labels   = {'A6', 'D6', 'D5', 'D4', 'D3', 'D2', 'D1'};
total_energy_V = sum(Vload.^2);
total_energy_I = sum(Iout.^2);

fprintf('\n===========================================================================\n');
fprintf('     COEFFICIENT ENERGY AND PERCENTAGE DISTRIBUTION (VOLTAGE & CURRENT)      \n');
fprintf('===========================================================================\n');
fprintf('%-6s | %-10s | %-14s | %-10s || %-14s | %-10s\n', 'Coeff', 'Length', 'V Energy', '% V Energy', 'I Energy', '% I Energy');
fprintf('---------------------------------------------------------------------------\n');
for i = 1:length(coeffs_V)
    e_V = sum(coeffs_V{i}.^2);
    p_V = (e_V / total_energy_V) * 100;
    
    e_I = sum(coeffs_I{i}.^2);
    p_I = (e_I / total_energy_I) * 100;
    
    fprintf('%-6s | %-10d | %-14.4e | %-9.4f%% || %-14.4e | %-9.4f%%\n', ...
        labels{i}, length(coeffs_V{i}), e_V, p_V, e_I, p_I);
end
fprintf('---------------------------------------------------------------------------\n');

% =========================================================
%      DYNAMIC LOAD IMPEDANCE PLOT (200 ms WINDOWS)
% =========================================================
% Time vector corresponding to the midpoint of each 200 ms window
t_win = ((0:num_windows-1) + 0.5) * Twin; 

figure('Name', 'Dynamic Load Impedance Tracking', 'NumberTitle', 'off');
plot(t_win, R, '-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Resistance R (Real)');
hold on;
plot(t_win, X, '-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Reactance X (Imaginary)');
hold off;

grid on;
xlabel('Time (s)');
ylabel('Impedance (\Omega)');
title('Dynamic Load Impedance Tracking Z_{load} (200 ms Windows)');
legend('Location', 'best');
xlim([0, T]);

% =========================================================
%      ADDED: ACTUAL LOAD IMPEDANCE AND RESULTS SUMMARY
% =========================================================
% --- EDIT TO MATCH YOUR SIMULINK CIRCUIT ---------------------------------
% One entry per parallel load branch (series R + L + C within each branch).
% Use 0 for L or C if the branch does not contain that element.
R_br  = [1        10   2         5       ];   % ohm
L_br  = [3.18e-3  0    31.8e-3   0       ];   % H
C_br  = [0        0    0         15.92e-3];   % F
t_on  = [0        2.5  5.0       7.5     ];   % time each branch connects (s)
f0     = 48;                                  % source frequency at t = 0 (Hz)
f_rate = 0.4;                                 % frequency ramp (Hz/s): 48 -> 52 Hz over 10 s
guard  = 0.1;                                 % windows within this time of a switch are treated as transient (s)
% -------------------------------------------------------------------------

% Actual impedance evaluated at the instantaneous source frequency
actual_Z = @(tt) 1 ./ sum(cell2mat(arrayfun(@(b) ...
    (tt(:) >= t_on(b)) ./ (R_br(b) + 1i*2*pi*(f0 + f_rate*tt(:))*L_br(b) + ...
    (C_br(b) > 0) ./ (1i*2*pi*(f0 + f_rate*tt(:))*max(C_br(b), eps))), ...
    1:numel(R_br), 'UniformOutput', false)), 2);

t_act = linspace(0, T, 4001)';
Z_act = actual_Z(t_act);

hold on;
plot(t_act, real(Z_act), '--', 'LineWidth', 1.2, 'Color', [0 0.447 0.741], 'DisplayName', 'Actual R');
plot(t_act, imag(Z_act), '--', 'LineWidth', 1.2, 'Color', [0.850 0.325 0.098], 'DisplayName', 'Actual X');
hold off;

% Actual values at each window midpoint, and which windows are steady state
Z_act_win = actual_Z(t_win(:));
R_act_win = real(Z_act_win);
X_act_win = imag(Z_act_win);
t_sw = t_on(t_on > 0);
steady = true(num_windows, 1);
for kk = 1:num_windows
    w_start = (kk-1) * Twin - guard;
    w_end   = kk * Twin + guard;
    if any(t_sw > w_start & t_sw < w_end)
        steady(kk) = false;
    end
end

dR = R(:) - R_act_win;
dX = X(:) - X_act_win;

fprintf('\n=====================================================================\n');
fprintf('                 IMPEDANCE TRACKING SUMMARY (steady state)            \n');
fprintf('=====================================================================\n');
fprintf('%-13s | %-8s %-12s | %-8s %-12s\n', 'Stage', 'R act', 'Max R err', 'X act', 'Max X err');
fprintf('---------------------------------------------------------------------\n');
pR = 100 * abs(dR) ./ abs(R_act_win);              % per-window % error in R
pX = 100 * abs(dX) ./ max(abs(X_act_win), eps);    % per-window % error in X
t_edges = [t_on, T];
for ss = 1:numel(t_on)
    in_stage = steady & t_win(:) >= t_edges(ss) & t_win(:) < t_edges(ss+1);
    if ~any(in_stage), continue; end
    Ra = mean(R_act_win(in_stage));
    Xa = mean(X_act_win(in_stage));
    fprintf('%4.1f - %4.1f s | %8.4f %10.3f%%  | %8.4f %10.3f%%\n', ...
        t_edges(ss), t_edges(ss+1), Ra, max(pR(in_stage)), Xa, max(pX(in_stage)));
end
fprintf('---------------------------------------------------------------------\n');
fprintf('RMSE  R: %.4f ohm   X: %.4f ohm   (%d of %d windows, excl. %.1f s around switches)\n', ...
    sqrt(mean(dR(steady).^2)), sqrt(mean(dX(steady).^2)), nnz(steady), num_windows, guard);
fprintf('Max |error|  R: %.4f ohm (%.3f%%)   X: %.4f ohm (%.3f%%)\n', ...
    max(abs(dR(steady))), max(pR(steady)), max(abs(dX(steady))), max(pX(steady)));
fprintf('Max |error| incl. switching windows  R: %.4f ohm   X: %.4f ohm\n', max(abs(dR)), max(abs(dX)));
fprintf('=====================================================================\n');
 