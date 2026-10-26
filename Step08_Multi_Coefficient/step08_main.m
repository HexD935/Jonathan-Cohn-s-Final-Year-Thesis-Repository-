Fs = 8192;
T  = 10;
N  = Fs * T;             
t  = (0:N-1)' / Fs;     

Vload  = out.simout.Data(:,1);
Iout   = out.simout.Data(:,2);

LoD = [ ...
   -0.07576571478927333, -0.02963552764599851,  0.49761866763201545,  0.80373875180591610, ...
    0.29785779560527736, -0.09921954357684722, -0.01260396726203783,  0.03222310060404270 ];
HiD = [ ...
   -0.03222310060404270, -0.01260396726203783,  0.09921954357684722,  0.29785779560527736, ...
   -0.80373875180591610,  0.49761866763201545,  0.02963552764599851, -0.07576571478927333 ];

% =========================================================================
%  ADAPTIVE 10-CYCLE SLIDING WINDOW SETUP (20 ms Update Step)
% =========================================================================
T_step = 0.020;                      % Fixed 20 ms update step
max_Tw = 10 / 48.0;                  % Max window duration (~208.33 ms at 48 Hz)

t_update = (0 : T_step : (T - max_Tw))'; 
num_updates = length(t_update);

f_search = 48.0 : 0.1 : 52.0;        % Ramping Fourier search band

% Preallocate storage for results across the 3 methods
f_est      = zeros(num_updates, 1);
Tw_vec     = zeros(num_updates, 1);

Z_A6       = zeros(num_updates, 1, 'like', 1i);
Z_A6D6D5   = zeros(num_updates, 1, 'like', 1i);
Z_TD       = zeros(num_updates, 1, 'like', 1i);

% =========================================================================
%  SLIDING WINDOW PROCESSING LOOP
% =========================================================================
[A6_Vf, D6_Vf, D5_Vf] = dwt_sym4_l6(Vload, LoD, HiD);
[A6_If, D6_If, D5_If] = dwt_sym4_l6(Iout,  LoD, HiD);
for k = 1:num_updates
    t_start = t_update(k);
    idx_start = round(t_start * Fs) + 1;
    
    % 1. Estimate fundamental frequency over nominal 200 ms interval
    idx_nom_end = min(N, idx_start + round(0.200 * Fs) - 1);
    V_search = Vload(idx_start:idx_nom_end);
    t_search = t(idx_start:idx_nom_end);
    
    E_f = zeros(length(f_search), 1);
    for fi = 1:length(f_search)
        f_cand = f_search(fi);
        X_f = sum(V_search .* exp(-1i * 2 * pi * f_cand * t_search));
        E_f(fi) = abs(X_f)^2;
    end
    [~, max_idx] = max(E_f);
    fk = f_search(max_idx);
    f_est(k) = fk;
    
    % 2. Calculate adaptive 10-cycle window duration (Tw = 10 / fk)
    Tw_k = 10 / fk;
    Tw_vec(k) = Tw_k;
    
    % 3. Define exact 10-cycle window sample bounds
    idx_end = min(N, idx_start + round(Tw_k * Fs) - 1);
    
    % 4. Construct window gate and apply to V, I, and complex reference q
    g = zeros(N, 1);
    g(idx_start:idx_end) = 1;
    
    q_k     = g .* exp(-1i * 2 * pi * fk * t);
    V_gated = g .* Vload;
    I_gated = g .* Iout;
    
    % 5. Execute 6-Level DWT (sym4) to extract A6, D6, and D5
    [A6_q, D6_q, D5_q] = dwt_sym4_l6(q_k,     LoD, HiD);
    
    % --- METHOD 1: A6-Only ---
    Vk_A6 = sum(A6_Vf .* A6_q);
    Ik_A6 = sum(A6_If .* A6_q);
    Z_A6(k) = Vk_A6 / Ik_A6;
    
    % --- METHOD 2: A6 + D6 + D5 ---
    Vk_A6D6D5 = sum(A6_Vf .* A6_q) + sum(D6_Vf .* D6_q) + sum(D5_Vf .* D5_q);
    Ik_A6D6D5 = sum(A6_If .* A6_q) + sum(D6_If .* D6_q) + sum(D5_If .* D5_q);
    Z_A6D6D5(k) = Vk_A6D6D5 / Ik_A6D6D5;
    
    % --- METHOD 3: Direct Time-Domain Reference ---
    Vk_TD = sum(V_gated .* q_k);
    Ik_TD = sum(I_gated .* q_k);
    Z_TD(k) = Vk_TD / Ik_TD;
end

R_A6     = real(Z_A6);     X_A6     = imag(Z_A6);
R_A6D6D5 = real(Z_A6D6D5); X_A6D6D5 = imag(Z_A6D6D5);
R_TD     = real(Z_TD);     X_TD     = imag(Z_TD);

% =========================================================================
%  GROUND TRUTH CALCULATIONS
% =========================================================================
t_mid  = t_update + (Tw_vec / 2);
f_true = 48 + (52 - 48) * (t_mid / T);

R1 = 1; L1 = 0.00318;            % Branch 1: 1 Ohm + 3.18 mH
R2 = 10;                         % Branch 2: 10 Ohm
R3 = 2; L3 = 0.0318;             % Branch 3: 2 Ohm + 31.8 mH
R4 = 5; C4 = 0.01592;            % Branch 4: 5 Ohm + 15.92 mF

Z_true = zeros(num_updates, 1, 'like', 1i);
for k = 1:num_updates
    tm = t_mid(k);
    fk = f_true(k);
    w  = 2 * pi * fk;
    
    Zb1 = R1 + 1i * w * L1;
    Zb2 = R2;
    Zb3 = R3 + 1i * w * L3;
    Zb4 = R4 - 1i / (w * C4);
    
    if tm < 2.5
        Y_total = 1 / Zb1;
    elseif tm < 5.0
        Y_total = 1 / Zb1 + 1 / Zb2;
    elseif tm < 7.5
        Y_total = 1 / Zb1 + 1 / Zb2 + 1 / Zb3;
    else
        Y_total = 1 / Zb1 + 1 / Zb2 + 1 / Zb3 + 1 / Zb4;
    end
    Z_true(k) = 1 / Y_total;
end
R_true = real(Z_true);
X_true = imag(Z_true);

% =========================================================================
%  QUANTITATIVE ERROR METRICS AND MASKING SWITCHING TRANSIENTS
% =========================================================================
% Define steady-state mask excluding 220 ms windows surrounding switching
switch_times = [2.5, 5.0, 7.5];
is_steady = true(num_updates, 1);
for st = switch_times
    is_steady = is_steady & ~(t_update >= (st - 0.20) & t_update <= (st + 0.22));
end

% Frequency metrics
f_error     = f_est - f_true;
f_rmse      = sqrt(mean(f_error.^2));
f_max_error = max(abs(f_error));

% Window duration metrics
Tw_min  = min(Tw_vec) * 1000; % ms
Tw_max  = max(Tw_vec) * 1000; % ms
Tw_mean = mean(Tw_vec) * 1000; % ms

% Impedance RMS Errors (Overall vs Steady-State)
calc_rmse = @(calc, true_val) sqrt(mean((calc - true_val).^2));

% Overall RMS Errors
rmse_R_A6_all     = calc_rmse(R_A6, R_true);
rmse_X_A6_all     = calc_rmse(X_A6, X_true);
rmse_R_Combo_all  = calc_rmse(R_A6D6D5, R_true);
rmse_X_Combo_all  = calc_rmse(X_A6D6D5, X_true);
rmse_R_TD_all     = calc_rmse(R_TD, R_true);
rmse_X_TD_all     = calc_rmse(X_TD, X_true);

% Steady-State RMS Errors (Excluding Switching Events)
rmse_R_A6_ss    = calc_rmse(R_A6(is_steady), R_true(is_steady));
rmse_X_A6_ss    = calc_rmse(X_A6(is_steady), X_true(is_steady));
rmse_R_Combo_ss = calc_rmse(R_A6D6D5(is_steady), R_true(is_steady));
rmse_X_Combo_ss = calc_rmse(X_A6D6D5(is_steady), X_true(is_steady));
rmse_R_TD_ss    = calc_rmse(R_TD(is_steady), R_true(is_steady));
rmse_X_TD_ss    = calc_rmse(X_TD(is_steady), X_true(is_steady));

% Print Statistical Summary
fprintf('\n===========================================================================\n');
fprintf('                  QUANTITATIVE PERFORMANCE SUMMARY                        \n');
fprintf('===========================================================================\n');
fprintf('Frequency Tracking RMSE    : %.4f Hz (Max Error: %.4f Hz)\n', f_rmse, f_max_error);
fprintf('Adaptive Window Duration   : Range = [%.2f ms, %.2f ms], Mean = %.2f ms\n', Tw_min, Tw_max, Tw_mean);
fprintf('---------------------------------------------------------------------------\n');
fprintf('%-25s | %-12s | %-12s | %-12s | %-12s\n', 'Method', 'R RMSE (All)', 'X RMSE (All)', 'R RMSE (SS)', 'X RMSE (SS)');
fprintf('---------------------------------------------------------------------------\n');
fprintf('%-25s | %-12.6f | %-12.6f | %-12.6f | %-12.6f\n', '1. A6 Only', rmse_R_A6_all, rmse_X_A6_all, rmse_R_A6_ss, rmse_X_A6_ss);
fprintf('%-25s | %-12.6f | %-12.6f | %-12.6f | %-12.6f\n', '2. A6 + D6 + D5', rmse_R_Combo_all, rmse_X_Combo_all, rmse_R_Combo_ss, rmse_X_Combo_ss);
fprintf('%-25s | %-12.6f | %-12.6f | %-12.6f | %-12.6f\n', '3. Direct Time-Domain', rmse_R_TD_all, rmse_X_TD_all, rmse_R_TD_ss, rmse_X_TD_ss);
fprintf('===========================================================================\n');

% =========================================================================
%  PLOTTING GRAPHICAL DELIVERABLES
% =========================================================================

% --- FIGURE 1: FREQUENCY & WINDOW DURATION ---
figure('Name', 'Frequency Tracking and Adaptive Window Duration', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(t_update, f_est, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5, 'DisplayName', 'Tracked Frequency f_k');
hold on;
plot(t_update, f_true, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True Ramp (48-52 Hz)');
hold off;
grid on; ylabel('Frequency (Hz)'); title('Estimated vs. Known Fundamental Frequency');
legend('Location', 'best'); xlim([0, T]);

subplot(2, 1, 2);
plot(t_update, Tw_vec * 1000, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5);
grid on; xlabel('Time (s)'); ylabel('Duration T_w (ms)');
title('Adaptive 10-Cycle Window Duration T_w = 10 / f_k'); xlim([0, T]);

% --- FIGURE 2: CALCULATED VS THEORETICAL IMPEDANCE (R & X) ---
figure('Name', 'Calculated vs. Theoretical Load Impedance (R & X)', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(t_update, R_A6D6D5, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', 'R (A6+D6+D5)');
hold on;
plot(t_update, R_true, 'k--', 'LineWidth', 2.0, 'DisplayName', 'R Theoretical');
hold off;
grid on; ylabel('Resistance (\Omega)'); title('Calculated vs. Theoretical Resistance R(t)');
legend('Location', 'best'); xlim([0, T]);

subplot(2, 1, 2);
plot(t_update, X_A6D6D5, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'DisplayName', 'X (A6+D6+D5)');
hold on;
plot(t_update, X_true, 'k--', 'LineWidth', 2.0, 'DisplayName', 'X Theoretical');
hold off;
grid on; xlabel('Time (s)'); ylabel('Reactance (\Omega)'); title('Calculated vs. Theoretical Reactance X(t)');
legend('Location', 'best'); xlim([0, T]);

% --- FIGURE 3: METHOD COMPARISON (A6-Only vs A6+D6+D5) ---
figure('Name', 'Method Comparison: A6-Only vs. A6+D6+D5', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(t_update, R_A6, 'r:', 'LineWidth', 1.5, 'DisplayName', 'R (A6 Only)');
hold on;
plot(t_update, R_A6D6D5, 'b-', 'LineWidth', 1.5, 'DisplayName', 'R (A6+D6+D5)');
plot(t_update, R_true, 'k--', 'LineWidth', 1.5, 'DisplayName', 'R Theoretical');
hold off;
grid on; ylabel('Resistance (\Omega)'); title('Resistance: A6-Only vs. A6+D6+D5');
legend('Location', 'best'); xlim([0, T]);

subplot(2, 1, 2);
plot(t_update, X_A6, 'r:', 'LineWidth', 1.5, 'DisplayName', 'X (A6 Only)');
hold on;
plot(t_update, X_A6D6D5, 'b-', 'LineWidth', 1.5, 'DisplayName', 'X (A6+D6+D5)');
plot(t_update, X_true, 'k--', 'LineWidth', 1.5, 'DisplayName', 'X Theoretical');
hold off;
grid on; xlabel('Time (s)'); ylabel('Reactance (\Omega)'); title('Reactance: A6-Only vs. A6+D6+D5');
legend('Location', 'best'); xlim([0, T]);

% --- FIGURE 4: METHOD COMPARISON (A6+D6+D5 vs Direct Time-Domain) ---
figure('Name', 'Method Comparison: A6+D6+D5 vs. Direct Time-Domain', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(t_update, R_A6D6D5, 'b-', 'LineWidth', 1.5, 'DisplayName', 'R (A6+D6+D5)');
hold on;
plot(t_update, R_TD, 'g--', 'LineWidth', 1.5, 'DisplayName', 'R (Direct Time-Domain)');
plot(t_update, R_true, 'k:', 'LineWidth', 1.5, 'DisplayName', 'R Theoretical');
hold off;
grid on; ylabel('Resistance (\Omega)'); title('Resistance: A6+D6+D5 vs. Direct Time-Domain');
legend('Location', 'best'); xlim([0, T]);

subplot(2, 1, 2);
plot(t_update, X_A6D6D5, 'b-', 'LineWidth', 1.5, 'DisplayName', 'X (A6+D6+D5)');
hold on;
plot(t_update, X_TD, 'g--', 'LineWidth', 1.5, 'DisplayName', 'X (Direct Time-Domain)');
plot(t_update, X_true, 'k:', 'LineWidth', 1.5, 'DisplayName', 'X Theoretical');
hold off;
grid on; xlabel('Time (s)'); ylabel('Reactance (\Omega)'); title('Reactance: A6+D6+D5 vs. Direct Time-Domain');
legend('Location', 'best'); xlim([0, T]);

% =========================================================================
%  LOCAL HELPER FUNCTION: 6-LEVEL SYM4 DWT DECONSTRUCTION
% =========================================================================
function [A6, D6, D5] = dwt_sym4_l6(sig, LoD, HiD)
    max_level = 6;
    A_curr = sig;
    D_coeff = cell(max_level, 1);
    
    for j = 1:max_level
        A_shifted = circshift(A_curr, -7);
        N_prev = length(A_shifted);
        N_curr = N_prev / 2;
        
        A_next = zeros(N_curr, 1, 'like', sig);
        D_next = zeros(N_curr, 1, 'like', sig);
        
        for n = 0:(N_curr - 1)
            sumA = 0; sumD = 0;
            for k = 0:7
                idx = mod(2*n + k, N_prev) + 1;
                sumA = sumA + LoD(k+1) * A_shifted(idx);
                sumD = sumD + HiD(k+1) * A_shifted(idx);
            end
            A_next(n+1) = sumA;
            D_next(n+1) = sumD;
        end
        A_curr = A_next;
        D_coeff{j} = D_next;
    end
    A6 = A_curr;
    D6 = D_coeff{6};
    D5 = D_coeff{5};
end

% =========================================================================
%  ADDED: RESULTS ANALYSIS - STEP 8 PART 1 (A6 ONLY, ADAPTIVE WINDOW)
%  Place this section BEFORE the "LOCAL HELPER FUNCTION" block: MATLAB does
%  not allow script code after a local function.
%  For a comparison with Step 7 (one window every 200 ms), set T_step = 0.200
%  at the top of the script before running.
% =========================================================================
% Same steady-state definition as Steps 6 and 7: a window is excluded if a
% switching instant lies within 'guard' seconds of it.
t_sw    = [2.5 5.0 7.5];          % switching instants (s)
guard   = 0.1;                    % s
t_edges = [0, t_sw, T];

% Method analysed. Part 1 = A6 only. (For Part 2, use R_A6D6D5 / X_A6D6D5.)
Re = R_A6(:);   Xe = X_A6(:);

% Earlier steps, for comparison: maximum steady-state error per stage (%)
step6_maxR = [4.245 3.428 6.665 4.996];   step6_maxX = [4.421 4.915 5.211 7.619];
step7_maxR = [3.808 3.101 6.585 3.326];   step7_maxX = [5.891 4.802 5.865 3.688];
step6_rmse = [0.0260 0.0269];             % Step 6 steady-state RMSE [R X] (ohm)

tw = t_mid(:);    ft = f_true(:);   fe = f_est(:);
Rt = R_true(:);   Xt = X_true(:);
nW = num_updates;

steady = true(nW, 1);
for kk = 1:nW
    w_start = t_update(kk) - guard;
    w_end   = t_update(kk) + Tw_vec(kk) + guard;
    if any(t_sw > w_start & t_sw < w_end)
        steady(kk) = false;
    end
end

% ---------------- Frequency tracking ----------------
df      = fe - ft;
pf      = 100 * abs(df) ./ ft;
sweepN  = round(0.200 * Fs);
on_grid = abs(fe - round(ft*10)/10) < 1e-9;

fprintf('\n=====================================================================\n');
fprintf('           STEP 8 PART 1: FREQUENCY TRACKING ANALYSIS                 \n');
fprintf('=====================================================================\n');
fprintf('Update step          : %.0f ms  (%d windows, %d sweeps)\n', T_step*1000, nW, nW);
fprintf('Window length        : %.2f - %.2f ms (mean %.2f ms)\n', ...
    min(Tw_vec)*1000, max(Tw_vec)*1000, mean(Tw_vec)*1000);
fprintf('Cost per sweep       : %d candidates x %d samples = %d complex multiplications\n', ...
    numel(f_search), sweepN, numel(f_search)*sweepN);
fprintf('Max |f error|        : %.3f Hz (%.4f %%)\n', max(abs(df)), max(pf));
fprintf('Mean |f error|       : %.4f Hz per window (%.4f %%)\n', mean(abs(df)), mean(pf));
fprintf('RMS f error          : %.3f Hz\n', sqrt(mean(df.^2)));
fprintf('On nearest 0.1 Hz    : %d of %d windows\n', nnz(on_grid), nW);

% ---------------- Load impedance (A6 only) ----------------
dR = Re - Rt;   dX = Xe - Xt;
pR = 100 * abs(dR) ./ abs(Rt);
pX = 100 * abs(dX) ./ max(abs(Xt), eps);

fprintf('\n=====================================================================\n');
fprintf('   STEP 8 PART 1: LOAD IMPEDANCE, A6 ONLY (steady state, max error %%)  \n');
fprintf('=====================================================================\n');
fprintf('%-13s | %-24s | %-24s\n', '', 'Max R error (%)', 'Max X error (%)');
fprintf('%-13s | %-7s %-7s %-8s | %-7s %-7s %-8s\n', 'Stage', 'Step 6', 'Step 7', 'Step 8', 'Step 6', 'Step 7', 'Step 8');
fprintf('---------------------------------------------------------------------\n');
for ss = 1:numel(t_edges)-1
    in_stage = steady & tw >= t_edges(ss) & tw < t_edges(ss+1);
    if ~any(in_stage), continue; end
    fprintf('%4.1f - %4.1f s | %7.3f %7.3f %8.3f | %7.3f %7.3f %8.3f\n', ...
        t_edges(ss), t_edges(ss+1), ...
        step6_maxR(ss), step7_maxR(ss), max(pR(in_stage)), ...
        step6_maxX(ss), step7_maxX(ss), max(pX(in_stage)));
end
fprintf('---------------------------------------------------------------------\n');
fprintf('RMSE (steady, %d of %d windows)  R: %.4f ohm (Step 6: %.4f)   X: %.4f ohm (Step 6: %.4f)\n', ...
    nnz(steady), nW, sqrt(mean(dR(steady).^2)), step6_rmse(1), sqrt(mean(dX(steady).^2)), step6_rmse(2));
fprintf('Max |error| steady state         R: %.4f ohm (%.3f%%)   X: %.4f ohm (%.3f%%)\n', ...
    max(abs(dR(steady))), max(pR(steady)), max(abs(dX(steady))), max(pX(steady)));
fprintf('Max |error| incl. switching      R: %.4f ohm             X: %.4f ohm\n', ...
    max(abs(dR)), max(abs(dX)));
cc = corrcoef(abs(df(steady)), pX(steady));
fprintf('Correlation |f error| vs X error : %.2f\n', cc(1,2));
fprintf('=====================================================================\n');

% =========================================================================
%  STEP 8 PART 1 FIGURE: A6-only impedance, tracked frequency, window length
%  Place BEFORE the "LOCAL HELPER FUNCTION" block (after the analysis section).
% =========================================================================
tx = t_mid;                                   % plot every window at its centre

fig = figure('Name', 'Step 8 Part 1: Adaptive Window (A6 Only)', 'NumberTitle', 'off', ...
             'Units', 'centimeters', 'Position', [2 2 18 20], 'Color', 'w');
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Top: load impedance (A6 only) vs actual ---
ax1 = nexttile(tl);
plot(ax1, tx, R_A6,  '-', 'Color', [0 0.4470 0.7410],      'LineWidth', 1.2, 'DisplayName', 'R calculated');
hold(ax1, 'on');
plot(ax1, tx, R_true, '--', 'Color', [0 0.20 0.50],     'LineWidth', 1.8, 'DisplayName', 'R actual');
plot(ax1, tx, X_A6,  '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2, 'DisplayName', 'X calculated');
plot(ax1, tx, X_true, '--', 'Color', [0.50 0.10 0.00],  'LineWidth', 1.8, 'DisplayName', 'X actual');
hold(ax1, 'off');
grid(ax1, 'on');
ylabel(ax1, 'Impedance (\Omega)');
title(ax1, 'Load impedance Z_{Load} (A6 only, adaptive 10-cycle window)');
legend(ax1, 'Location', 'eastoutside');

% --- Middle: tracked vs actual frequency ---
ax2 = nexttile(tl);
plot(ax2, tx, f_est, '-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.2, 'DisplayName', 'Tracked');
hold(ax2, 'on');
plot(ax2, tx, f_true, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Actual');
hold(ax2, 'off');
grid(ax2, 'on');
ylabel(ax2, 'Frequency (Hz)');
title(ax2, 'Source frequency');
legend(ax2, 'Location', 'eastoutside');
ylim(ax2, [47.5 52.5]);

% --- Bottom: adaptive window length ---
ax3 = nexttile(tl);
stairs(ax3, tx, Tw_vec * 1000, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5, 'DisplayName', 'T_w = 10 / f_k');
grid(ax3, 'on');
xlabel(ax3, 'Time (s)');
ylabel(ax3, 'Window length (ms)');
title(ax3, 'Adaptive window length');
legend(ax3, 'Location', 'eastoutside');

linkaxes([ax1 ax2 ax3], 'x');
xlim(ax1, [0 T]);

exportgraphics(fig, 'Ass8_Part1.png', 'Resolution', 300);

% ---- Part 2: max steady-state error per stage for all three methods ----
methods = {R_A6, X_A6; R_A6D6D5, X_A6D6D5; R_TD, X_TD};
fprintf('\n%-13s | %-26s | %-26s\n', 'Stage', 'Max R error (%)', 'Max X error (%)');
fprintf('%-13s | %8s %8s %8s | %8s %8s %8s\n', '', 'A6', 'A6D6D5', 'Direct', 'A6', 'A6D6D5', 'Direct');
for ss = 1:numel(t_edges)-1
    in_stage = steady & tw >= t_edges(ss) & tw < t_edges(ss+1);
    mR = zeros(1,3);  mX = zeros(1,3);
    for mm = 1:3
        pRm = 100 * abs(methods{mm,1}(:) - Rt) ./ abs(Rt);
        pXm = 100 * abs(methods{mm,2}(:) - Xt) ./ max(abs(Xt), eps);
        mR(mm) = max(pRm(in_stage));  mX(mm) = max(pXm(in_stage));
    end
    fprintf('%4.1f - %4.1f s | %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f\n', ...
        t_edges(ss), t_edges(ss+1), mR, mX);
end