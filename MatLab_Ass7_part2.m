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
%  DYNAMIC 200ms GATING & ADAPTIVE FOURIER FREQUENCY SEARCH
% ============================================================
Twin = 0.2;                         % 200 ms interval
num_windows = round(T / Twin);      % 50 windows total

f_search = 48.0 : 0.1 : 52.0;       % Search band: 48.0 to 52.0 Hz in 0.1 Hz steps
freq_update_interval = 1;           % Run search every 2nd window (change to 3 for 3rd)

q = cell(num_windows, 1);
f_tracked = zeros(num_windows, 1);  % Track detected frequency per window
f_current = 50.0;                   % Initial default frequency

for k = 1:num_windows
    idx_start = round((k-1) * Fs * Twin) + 1;
    idx_end   = round(k * Fs * Twin);
    
    % Update estimated fundamental frequency periodically using Voltage signal
    if mod(k-1, freq_update_interval) == 0
        V_win = Vload(idx_start:idx_end);
        t_win = t(idx_start:idx_end);
        
        E_f = zeros(length(f_search), 1);
        for fi = 1:length(f_search)
            f_cand = f_search(fi);
            % Ramping Fourier transform integral over voltage window
            X_f = sum(V_win .* exp(-1i * 2 * pi * f_cand * t_win));
            E_f(fi) = abs(X_f)^2;   % Signal energy at candidate frequency
        end
        
        % Extract frequency corresponding to maximum energy peak
        [~, max_idx] = max(E_f);
        f_current = f_search(max_idx);
    end
    
    f_tracked(k) = f_current;
    
    % Construct window gate and update reference signal with tracked frequency
    g = zeros(N, 1);
    g(idx_start:idx_end) = 1;
    q{k} = g .* exp(-1i * 2 * pi * f_current * t);
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

% Compute Phasors and Impedance for each 200 ms window
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
%      GROUND TRUTH (REFERENCE) CALCULATIONS
% =========================================================
t_win = ((0:num_windows-1) + 0.5) * Twin; % Midpoint time for 200 ms windows

% 1. Ground Truth Frequency Ramp (48 Hz to 52 Hz over T = 10s)
f_true = 48 + (52 - 48) * (t_win / T);

% 2. Component Values for Dynamic Load Impedances
R1 = 1; L1 = 0.00318;            % Branch 1: 1 Ohm + 0.00318 H
R2 = 10;                         % Branch 2: 10 Ohm
R3 = 2; L3 = 0.0318;             % Branch 3: 2 Ohm + 0.0318 H
R4 = 5; C4 = 0.01592;            % Branch 4: 5 Ohm + 0.01592 F

Z_true = zeros(num_windows, 1, 'like', 1i);

for k = 1:num_windows
    tk = t_win(k);
    fk = f_true(k);
    w  = 2 * pi * fk;            % Frequency-dependent angular frequency
    
    % Dynamic Branch Impedances
    Zb1 = R1 + 1i * w * L1;
    Zb2 = R2;
    Zb3 = R3 + 1i * w * L3;
    Zb4 = R4 - 1i / (w * C4);
    
    % Parallel Combinations across Stages
    if tk < 2.5
        Y_total = 1 / Zb1;
    elseif tk < 5.0
        Y_total = 1 / Zb1 + 1 / Zb2;
    elseif tk < 7.5
        Y_total = 1 / Zb1 + 1 / Zb2 + 1 / Zb3;
    else
        Y_total = 1 / Zb1 + 1 / Zb2 + 1 / Zb3 + 1 / Zb4;
    end
    
    Z_true(k) = 1 / Y_total;
end

R_true = real(Z_true);
X_true = imag(Z_true);

% =========================================================
%      DYNAMIC LOAD IMPEDANCE AND FREQUENCY PLOTS
% =========================================================
figure('Name', 'Dynamic Load Impedance and Frequency Tracking with Ground Truth', 'NumberTitle', 'off');

% Subplot 1: Dynamic Load Impedance vs. Ground Truth
subplot(2, 1, 1);
plot(t_win, R, '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'R Calculated');
hold on;
plot(t_win, R_true, '--', 'Color', [0 0.20 0.50], 'LineWidth', 2.0, 'DisplayName', 'R True');

plot(t_win, X, '-s', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'X Calculated');
plot(t_win, X_true, '--', 'Color', [0.50 0.10 0.00], 'LineWidth', 2.0, 'DisplayName', 'X True');
hold off;

grid on;
ylabel('Impedance (\Omega)');
title('Dynamic Load Impedance Tracking Z_{load} vs. Ground Truth (200 ms Windows)');
legend('Location', 'best');
xlim([0, T]);

% Subplot 2: Tracked Fundamental Frequency vs. Ground Truth Ramp
subplot(2, 1, 2);
plot(t_win, f_tracked, '-d', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Frequency Tracked');
hold on;
plot(t_win, f_true, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Frequency True (Ramp 48-52 Hz)');
hold off;

grid on;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Fundamental Frequency Search vs. Ground Truth Ramp');
legend('Location', 'best');
ylim([47.0, 53.0]);
xlim([0, T]);

% =========================================================
%      ADDED: RESULTS ANALYSIS (FREQUENCY AND IMPEDANCE)
% =========================================================
% Same steady-state definition as the Step 6 analysis, so results compare:
% windows within 'guard' seconds of a switch are excluded.
t_sw  = [2.5 5.0 7.5];            % switching instants (s), as in the ground truth above
guard = 0.1;                      % s
t_edges = [0, t_sw, T];

% Step 6 results (fixed 50 Hz reference), for comparison only
step6_maxR = [4.777 4.197 4.849 4.234];   % max R error per stage (%)
step6_maxX = [8.661 6.127 5.791 6.567];   % max X error per stage (%)
step6_rmse = [0.0280 0.0323];             % steady-state RMSE [R X] (ohm)

tw  = t_win(:);
ft  = f_true(:);
fe  = f_tracked(:);
Rt  = R_true(:);  Xt = X_true(:);
Re  = R(:);       Xe = X(:);

steady = true(num_windows, 1);
for kk = 1:num_windows
    w_start = (kk-1) * Twin - guard;
    w_end   = kk * Twin + guard;
    if any(t_sw > w_start & t_sw < w_end)
        steady(kk) = false;
    end
end

% ---------------- Frequency tracking ----------------
df       = fe - ft;
n_sweeps = numel(1:freq_update_interval:num_windows);
win_len  = round(Fs * Twin);
on_grid  = abs(fe - round(ft*10)/10) < 1e-9;   % estimate equals nearest 0.1 Hz grid point

fprintf('\n=====================================================================\n');
fprintf('                    FREQUENCY TRACKING ANALYSIS                       \n');
fprintf('=====================================================================\n');
fprintf('Update interval      : every %d window(s)  (%d sweeps over %d windows)\n', ...
    freq_update_interval, n_sweeps, num_windows);
fprintf('Cost per sweep       : %d candidates x %d samples = %d complex multiplications\n', ...
    numel(f_search), win_len, numel(f_search)*win_len);
fprintf('Max |f error|        : %.3f Hz (all windows)   %.3f Hz (steady state)\n', ...
    max(abs(df)), max(abs(df(steady))));
fprintf('RMS f error          : %.3f Hz (all windows)   %.3f Hz (steady state)\n', ...
    sqrt(mean(df.^2)), sqrt(mean(df(steady).^2)));
fprintf('Mean f error (bias)  : %+.3f Hz\n', mean(df));
pf = 100 * abs(df) ./ ft;                       % per-window % frequency error
fprintf('Max %% f error        : %.4f %% (over the 10 s)\n', max(pf));
fprintf('Mean |f error|       : %.4f Hz per window  (%.4f %%)\n', mean(abs(df)), mean(pf));
fprintf('On nearest 0.1 Hz    : %d of %d windows\n', nnz(on_grid), num_windows);
[~, kf] = max(abs(df));
fprintf('Worst window         : %.1f s (tracked %.1f Hz, actual %.2f Hz)\n', tw(kf), fe(kf), ft(kf));

% ---------------- Load impedance ----------------
dR = Re - Rt;   dX = Xe - Xt;
pR = 100 * abs(dR) ./ abs(Rt);
pX = 100 * abs(dX) ./ max(abs(Xt), eps);

fprintf('\n=====================================================================\n');
fprintf('             LOAD IMPEDANCE ANALYSIS (steady state, max error)        \n');
fprintf('=====================================================================\n');
fprintf('%-13s | %-7s %-7s | %-16s | %-16s\n', '', 'R act', 'X act', 'Max R err (%)', 'Max X err (%)');
fprintf('%-13s | %-7s %-7s | %-7s %-8s | %-7s %-8s\n', 'Stage', '', '', 'Step 6', 'Step 7', 'Step 6', 'Step 7');
fprintf('---------------------------------------------------------------------\n');
for ss = 1:numel(t_edges)-1
    in_stage = steady & tw >= t_edges(ss) & tw < t_edges(ss+1);
    if ~any(in_stage), continue; end
    fprintf('%4.1f - %4.1f s | %7.4f %7.4f | %7.3f %8.3f | %7.3f %8.3f\n', ...
        t_edges(ss), t_edges(ss+1), mean(Rt(in_stage)), mean(Xt(in_stage)), ...
        step6_maxR(ss), max(pR(in_stage)), step6_maxX(ss), max(pX(in_stage)));
end
fprintf('---------------------------------------------------------------------\n');
rmseR = sqrt(mean(dR(steady).^2));  rmseX = sqrt(mean(dX(steady).^2));
fprintf('RMSE (steady, %d of %d windows)  R: %.4f ohm (Step 6: %.4f)   X: %.4f ohm (Step 6: %.4f)\n', ...
    nnz(steady), num_windows, rmseR, step6_rmse(1), rmseX, step6_rmse(2));
fprintf('Max |error| steady state         R: %.4f ohm (%.3f%%)   X: %.4f ohm (%.3f%%)\n', ...
    max(abs(dR(steady))), max(pR(steady)), max(abs(dX(steady))), max(pX(steady)));
fprintf('Max |error| incl. switching      R: %.4f ohm             X: %.4f ohm\n', ...
    max(abs(dR)), max(abs(dX)));
[~, kx] = max(pX .* steady);
fprintf('Worst steady X window            : %.1f s (%.3f%%, f error %+.3f Hz)\n', tw(kx), pX(kx), df(kx));

% Does the remaining X error follow the frequency error?
cc = corrcoef(abs(df(steady)), pX(steady));
fprintf('Correlation |f error| vs X error : %.2f  (near 0 = X error not driven by frequency)\n', cc(1,2));
fprintf('=====================================================================\n');