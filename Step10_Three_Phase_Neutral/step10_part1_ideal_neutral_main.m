% =========================================================================
%  MULTI-COEFFICIENT (SYM4 DWT) IMPEDANCE ESTIMATION - THREE PHASES
%    1. Z_load of every phase from the sliding-window DWT coefficients
%    2. Switch detection and Thevenin impedance Zth of every phase
%  The Simulink output out.simout.Data has the columns  V1 I1 V2 I2 V3 I3.
%  The three phases share one window and one frequency estimate per update;
%  everything else (Z_load, switch detection, Zth) is done phase by phase.
% =========================================================================
Fs = 8192;
T  = 10;
N  = Fs * T;
t  = (0:N-1)' / Fs;

n_ph = 3;                                 % number of phases
if size(out.simout.Data, 2) < 2 * n_ph
    error('out.simout.Data must have the columns V1 I1 V2 I2 V3 I3 (found %d columns).', size(out.simout.Data, 2));
end
Vload = out.simout.Data(:, [1 3 5]);      % V1 V2 V3   (N x 3)
Iout  = out.simout.Data(:, [2 4 6]);      % I1 I2 I3   (N x 3)

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

% Preallocate storage for multi-coefficient calculations
f_est   = zeros(num_updates, 1);                  % one frequency estimate per window, shared by the three phases
Tw_vec  = zeros(num_updates, 1);
Nw_vec  = zeros(num_updates, 1);                  % number of samples in every window
tc_vec  = zeros(num_updates, 1);                  % exact centre time of every window (s)
Z_calc  = zeros(num_updates, n_ph, 'like', 1i);   % Z_load of every phase
Ik_vec  = zeros(num_updates, n_ph, 'like', 1i);   % current coefficient sum of every phase

% =========================================================================
%  CHANGED: DECOMPOSE V AND I OF EVERY PHASE ONCE OVER THE FULL RECORD
%  (not gated per window). Only the complex reference q is gated inside the
%  loop. Every window's phasor then matches the direct time-domain result,
%  and the loop needs one DWT per window instead of seven.
% =========================================================================
A6_Vf = zeros(N/2^6, n_ph);  D6_Vf = zeros(N/2^6, n_ph);  D5_Vf = zeros(N/2^5, n_ph);
A6_If = zeros(N/2^6, n_ph);  D6_If = zeros(N/2^6, n_ph);  D5_If = zeros(N/2^5, n_ph);
for ph = 1:n_ph
    [A6_Vf(:, ph), D6_Vf(:, ph), D5_Vf(:, ph)] = dwt_sym4_l6(Vload(:, ph), LoD, HiD);
    [A6_If(:, ph), D6_If(:, ph), D5_If(:, ph)] = dwt_sym4_l6(Iout(:, ph),  LoD, HiD);
end

% =========================================================================
%  SLIDING WINDOW PROCESSING LOOP
% =========================================================================
for k = 1:num_updates
    t_start = t_update(k);
    idx_start = round(t_start * Fs) + 1;

    % 1. Estimate fundamental frequency over nominal 200 ms interval
    %    (spectral energy of the three phase voltages added up)
    idx_nom_end = min(N, idx_start + round(0.200 * Fs) - 1);
    V_search = Vload(idx_start:idx_nom_end, :);
    t_search = t(idx_start:idx_nom_end);

    E_f = zeros(length(f_search), 1);
    for fi = 1:length(f_search)
        f_cand = f_search(fi);
        X_f = sum(V_search .* exp(-1i * 2 * pi * f_cand * t_search), 1);   % 1 x 3, one value per phase
        E_f(fi) = sum(abs(X_f).^2);
    end
    [~, max_idx] = max(E_f);
    fk = f_search(max_idx);
    f_est(k) = fk;

    % 2. Calculate adaptive 10-cycle window duration (Tw = 10 / fk)
    Tw_k = 10 / fk;
    Tw_vec(k) = Tw_k;

    % 3. Define exact 10-cycle window sample bounds
    idx_end = min(N, idx_start + round(Tw_k * Fs) - 1);

    % 4. Window gate and complex reference q (the same for the three phases); only q is gated
    g = zeros(N, 1);
    g(idx_start:idx_end) = 1;
    q_k = g .* exp(-1i * 2 * pi * fk * t);
    [A6_q, D6_q, D5_q] = dwt_sym4_l6(q_k, LoD, HiD);

    % 5. A6 + D6 + D5 coefficient sums of every phase
    %    (CHANGED: uses the full-record V and I coefficients from before the loop; no per-window DWT of V and I)
    for ph = 1:n_ph
        Vk_calc = sum(A6_Vf(:, ph) .* A6_q) + sum(D6_Vf(:, ph) .* D6_q) + sum(D5_Vf(:, ph) .* D5_q);
        Ik_calc = sum(A6_If(:, ph) .* A6_q) + sum(D6_If(:, ph) .* D6_q) + sum(D5_If(:, ph) .* D5_q);
        Z_calc(k, ph) = Vk_calc / Ik_calc;
        Ik_vec(k, ph) = Ik_calc;                     % kept for the Thevenin calculation
    end
    Nw_vec(k) = idx_end - idx_start + 1;             % samples in this window
    tc_vec(k) = (idx_start + idx_end - 2) / (2*Fs);  % centre time of this window (s)
end

R_calc = real(Z_calc);                               % num_updates x n_ph
X_calc = imag(Z_calc);

% =========================================================================
%  GROUND TRUTH CALCULATIONS
% =========================================================================
t_mid  = t_update + (Tw_vec / 2);
f_true = 48 + (52 - 48) * (t_mid / T);

% Load of every phase (ground truth, used only for the comparisons): four branches, each a series R-L-C, added in
% parallel at the switching instants. Columns = [phase 1, phase 2, phase 3]; L = 0: no inductor; C = Inf: no capacitor.
% Switching instants of every phase (one row per phase): branches 2, 3 and 4 are added at these times.
switch_times = [2.5  5.0  7.5;                    % phase 1
                2.0  6.0  8.0;                    % phase 2
                2.0  3.0  7.0];                   % phase 3
%        phase 1    phase 2   phase 3
R1 = [   1          0.8       0.6      ];   L1 = [ 0.00318   0         0.00191  ];   C1 = [ Inf       Inf       Inf      ];
R2 = [   10         0.2       1.5      ];   L2 = [ 0         0.00318   0        ];   C2 = [ Inf       Inf       Inf      ];
R3 = [   2          0.5       0.1      ];   L3 = [ 0.0318    0         0        ];   C3 = [ Inf       0.006366  0.003979 ];
R4 = [   5          2.0       0.3      ];   L4 = [ 0         0         0.001273 ];   C4 = [ 0.01592   Inf       Inf      ];
Z_true = zeros(num_updates, n_ph, 'like', 1i);

for ph = 1:n_ph
    for k = 1:num_updates
        w  = 2 * pi * f_true(k);
        Zb = [branch_z(R1(ph), L1(ph), C1(ph), w), ...
              branch_z(R2(ph), L2(ph), C2(ph), w), ...
              branch_z(R3(ph), L3(ph), C3(ph), w), ...
              branch_z(R4(ph), L4(ph), C4(ph), w)];
        n_on = 1 + sum(t_mid(k) >= switch_times(ph, :));   % branches connected at this time
        Z_true(k, ph) = 1 / sum(1 ./ Zb(1:n_on));    % parallel combination
    end
end

R_true = real(Z_true);
X_true = imag(Z_true);

% =========================================================================
%  QUANTITATIVE ERROR METRICS (per phase)
% =========================================================================
is_steady = true(num_updates, n_ph);              % steady-state windows of every phase (its own switching times)
for ph = 1:n_ph
    for st = switch_times(ph, :)
        is_steady(:, ph) = is_steady(:, ph) & ~(t_update >= (st - 0.20) & t_update <= (st + 0.22));
    end
end

calc_rmse = @(calc, true_val) sqrt(mean((calc - true_val).^2));
rmse_R_all = zeros(1, n_ph);  rmse_X_all = zeros(1, n_ph);
rmse_R_ss  = zeros(1, n_ph);  rmse_X_ss  = zeros(1, n_ph);
for ph = 1:n_ph
    ss = is_steady(:, ph);
    rmse_R_all(ph) = calc_rmse(R_calc(:, ph), R_true(:, ph));
    rmse_X_all(ph) = calc_rmse(X_calc(:, ph), X_true(:, ph));
    rmse_R_ss(ph)  = calc_rmse(R_calc(ss, ph), R_true(ss, ph));
    rmse_X_ss(ph)  = calc_rmse(X_calc(ss, ph), X_true(ss, ph));
end

fprintf('\n===========================================================================\n');
fprintf('         MULTI-COEFFICIENT METHOD PERFORMANCE SUMMARY (3 PHASES)            \n');
fprintf('===========================================================================\n');
for ph = 1:n_ph
    fprintf('Phase %d | Resistance RMS Error | Overall: %.6f Ohm | Steady-State: %.6f Ohm\n', ph, rmse_R_all(ph), rmse_R_ss(ph));
    fprintf('        | Reactance RMS Error  | Overall: %.6f Ohm | Steady-State: %.6f Ohm\n', rmse_X_all(ph), rmse_X_ss(ph));
end
fprintf('===========================================================================\n');

% =========================================================================
%  PLOTTING: ACTUAL VS CALCULATED Z_LOAD OF EACH PHASE (MULTI-COEFFICIENT METHOD)
% =========================================================================
figure('Name', 'Actual vs Calculated Impedance (Multi-Coefficient), 3 phases', 'NumberTitle', 'off', ...
       'Position', [60 60 1150 820]);
for ph = 1:n_ph
    % Resistance
    subplot(n_ph, 2, 2*ph - 1);
    plot(t_update, R_true(:, ph), 'k--', 'LineWidth', 2.0, 'DisplayName', 'Actual R');
    hold on;
    plot(t_update, R_calc(:, ph), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Calculated R (Multi-Coeff)');
    hold off;
    grid on;
    ylabel('Resistance (\Omega)');
    title(sprintf('Phase %d: Resistance, actual vs. calculated', ph));
    if ph == 1, legend('Location', 'best', 'FontSize', 8); end
    if ph == n_ph, xlabel('Time (s)'); end
    xlim([0, T]);

    % Reactance
    subplot(n_ph, 2, 2*ph);
    plot(t_update, X_true(:, ph), 'k--', 'LineWidth', 2.0, 'DisplayName', 'Actual X');
    hold on;
    plot(t_update, X_calc(:, ph), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Calculated X (Multi-Coeff)');
    hold off;
    grid on;
    ylabel('Reactance (\Omega)');
    title(sprintf('Phase %d: Reactance, actual vs. calculated', ph));
    if ph == 1, legend('Location', 'best', 'FontSize', 8); end
    if ph == n_ph, xlabel('Time (s)'); end
    xlim([0, T]);
end

% =========================================================================
%  SWITCH DETECTION + THEVENIN IMPEDANCE OF EACH PHASE (DWT-coefficient domain)
%
%      P   = I_1 / I_2                    1 = state after the transient
%      Zth = (Z_L2 - P*Z_L1) / (P - 1)    2 = state before the transient
%
%  Every phase is processed on its own by thevenin_from_load_changes() (bottom of
%  this file): its switches are detected from its own Z_calc (switch_times is NOT used
%  there) and only its coefficient-domain results (Ik_calc, Z_calc) enter the formula -
%  no waveform is rebuilt in the time domain.
% =========================================================================
% --- Settings (the same for the three phases, except where marked per phase) -------
% A switch is spread over one window length (~10 updates), so the change per update is only
% a fraction (~1/10) of the full load step. Set the thresholds above the steady-state
% update-to-update noise of Z_calc and below ~1/10 of the smallest load step to catch.
% Per phase [phase 1, phase 2, phase 3]: above that phase's steady-state noise, below ~1/10 of its smallest load step.
dR_thresh = [0.008 0.004 0.004];   % [Ohm per update] change of the REAL part between consecutive windows that flags a switch
dX_thresh = [0.008 0.004 0.004];   % [Ohm per update] change of the IMAGINARY part that flags a switch
n_guard   = 2;       % windows before the first flagged one that are discarded as well
t_settle  = 0.22;    % [s] windows start counting as "post-transient" this long after the switch has passed through
n_avg     = 5;       % windows averaged on each side of a switch (5 x 20 ms = 0.1 s)
% Source-phase smoother: frequency wander slower than f_smooth is followed. Higher follows faster
% wander but lets more measurement noise into the phase that is bridged across each switch; lower is smoother.
% 0.75 suits the linear ramp (on it the phase error of P is about 2 mrad, against about 9 mrad at 1.5 Hz in phases 2
% and 3); raise it to ~1.5 for random frequency wander - the printed fit residual tells when it is too stiff.
% (CHANGED: those figures were measured with the source phase tracked from I. It is now tracked from V (section 3 of
% thevenin_from_load_changes), which makes the phase error of P smaller; 0.75 Hz is still the better setting for the ramp.)
f_smooth  = 0.75;    % [Hz]

% Which part of Zth can a switch measure? An error in P reaches Zth through err_gain = dZ*P/(P-1)^2 (dZ = ZL1 - ZL2).
% The dominant error of P is a few mrad in its PHASE (bridging the source phase across the switch); that lands in
% the part of Zth the load change moved (R for a resistive change, X for a reactive one, with loads near 45 deg),
% so such a switch measures the OTHER part well. Each switch therefore updates Re(Zth) or Im(Zth) only if the
% predicted error of that part is small; parts measured by several switches are averaged with weight 1/error^2.
sig_phase = 3e-3;    % [rad] expected error of arg(P)
sig_amp   = 1e-3;    % [-]   expected relative error of |P|
tol_rel   = 0.05;    % accuracy required for a Zth change: a switch may update R_th (X_th) only if its predicted error
                     % is below this fraction of |Zth| (was 0.10; smaller = stricter, 0.03 or less rejects most switches)

% ACTUAL source (Thevenin) impedance of every phase [phase 1, phase 2, phase 3], Zth = Rth + j*2*pi*f*Lth.
% Phase 1: 0.1 + 0.2i. Phases 2 and 3: 0.03 Ohm + 0.000127 H and 0.05 Ohm + 0.000191 H (0.04 and 0.06 Ohm of reactance at 50 Hz).
% They are the dashed 'Actual' lines and the error figures, and the START VALUE of every Zth graph (below).
Rth = [0.1,  0.03,      0.05    ];          % [Ohm]
Lth = [0.2/(2*pi*50),  0.000127,  0.000191];   % [H]   (phase 1: 0.2 Ohm of reactance at 50 Hz)
f_nom     = 50;                                     % [Hz] nominal frequency: start value, and X <-> L conversion of the start value
Zth_start = Rth + 1i * 2 * pi * f_nom * Lth;        % [Ohm] where each Zth graph starts: the actual value at f_nom;
                                                    % a switch replaces it only when it measures that part accurately enough
Zth_true = Rth + 1i * 2 * pi * (f_true .* Lth);   % num_updates x n_ph

cfg = struct('dR_thresh', dR_thresh(1), 'dX_thresh', dX_thresh(1), 'n_guard', n_guard, 't_settle', t_settle, ...
             'n_avg', n_avg, 'f_smooth', f_smooth, 'sig_phase', sig_phase, 'sig_amp', sig_amp, ...
             'tol_rel', tol_rel, 'f_nom', f_nom, 'T', T, 'T_step', T_step, 'max_Tw', max_Tw);

res = cell(n_ph, 1);
for ph = 1:n_ph
    cfg.dR_thresh = dR_thresh(ph);                    % detection thresholds of this phase
    cfg.dX_thresh = dX_thresh(ph);
    res{ph} = thevenin_from_load_changes(Z_calc(:, ph), Ik_vec(:, ph), Nw_vec, tc_vec, f_est, Tw_vec, ...
                                         t_update, Zth_true(:, ph), Zth_start(ph), cfg);
end

% --- Report -------------------------------------------------------------------
fprintf('\n===========================================================================\n');
fprintf('  SWITCH DETECTION + THEVENIN IMPEDANCE OF EACH PHASE (P = I1/I2, DWT domain) \n');
fprintf('===========================================================================\n');
fprintf('A switch changes R_th / X_th only if its predicted error is below %.0f%% of |Zth| (the [limit] printed below)\n', 100 * tol_rel);
used_txt = {'not used', 'used'};
for ph = 1:n_ph
    r = res{ph};
    fprintf('\n--- Phase %d ----------------------------------------------------------------\n', ph);
    fprintf('%d event(s) detected: |dR| > %.4g Ohm or |dX| > %.4g Ohm between consecutive windows\n', ...
            r.n_ev, dR_thresh(ph), dX_thresh(ph));
    fprintf('Source-phase spline (phase tracked from V): f_smooth = %.2g Hz, fit residual %.2f mrad rms', f_smooth, r.fit_rms_mrad);   % CHANGED
    if r.fit_rms_mrad > 5, fprintf(' (above ~5 mrad: try a higher f_smooth)'); end
    fprintf('\n');
    fprintf('The Zth graph starts from the start value %.4f %+.4fj Ohm\n', real(Zth_start(ph)), imag(Zth_start(ph)));
    for e = 1:r.n_ev
        [dt_min, j] = min(abs(switch_times(ph, :) - r.t_event(e)));   % known times only used to show the detection error
        if dt_min < 0.1
            fprintf('Event %d at %.3f s (actual %.1f s)', e, r.t_event(e), switch_times(ph, j));
        else
            fprintf('Event %d at %.3f s (no known switch nearby)', e, r.t_event(e));
        end
        if isnan(r.Zth_event(e))
            fprintf(' | not used (no clean state before/after, or no net load change)\n');
            continue;
        end
        fprintf(' | |P-1| = %.3f\n', abs(r.P_event(e) - 1));
        fprintf('    R_th = %.4f +- %.4f [limit %.4f] (%s)    X_th = %.4f +- %.4f [limit %.4f] (%s)\n', ...
                real(r.Zth_event(e)), r.sig_R_event(e), tol_rel * abs(r.Zth_event(e)), used_txt{r.use_R(e) + 1}, ...
                imag(r.Zth_event(e)), r.sig_X_event(e), tol_rel * abs(r.Zth_event(e)), used_txt{r.use_X(e) + 1});
        fprintf('    L_th = %.4e H  (X_th at %.2f Hz%s)\n', r.L_event(e), r.f_event(e), ...
                sprintf('; actual %.4e H', Lth(ph)));
        if ~isnan(Rth(ph)) && ~isnan(Lth(ph))
            fprintf('    true  = %.4f %+.4fj | actual error: R %+.4f, X %+.4f Ohm\n', ...
                    real(r.Zth_true_event(e)), imag(r.Zth_true_event(e)), ...
                    real(r.Zth_event(e) - r.Zth_true_event(e)), imag(r.Zth_event(e) - r.Zth_true_event(e)));
        end
    end
    if any(r.use_R), txt_R = sprintf('%d switch(es) +- %.4f', sum(r.use_R), r.sig_R_final); else, txt_R = 'start value'; end
    if any(r.use_X), txt_X = sprintf('%d switch(es) +- %.4f', sum(r.use_X), r.sig_X_final); else, txt_X = 'start value'; end
    fprintf('Phase %d Zth = %.4f %+.4fj Ohm   (R_th: %s | X_th: %s)\n', ph, real(r.Zth_final), imag(r.Zth_final), txt_R, txt_X);
    if ~isnan(Rth(ph)) && ~isnan(Lth(ph))
        fprintf('Actual   Zth = %.4f %+.4fj Ohm   (error: R %+.4f, X %+.4f)\n', real(r.Zth_true_final), imag(r.Zth_true_final), ...
                real(r.Zth_final - r.Zth_true_final), imag(r.Zth_final - r.Zth_true_final));
    end
end

fprintf('\n===========================================================================\n');
fprintf('  SUMMARY OF THE THREE PHASES\n');
fprintf('===========================================================================\n');
fprintf('Phase | Z_load after the last switch [Ohm]: estimate / actual | Zth [Ohm]: estimate / actual\n');
for ph = 1:n_ph
    r   = res{ph};
    m   = r.is_clean & r.state == r.n_ev + 1;            % clean windows of the last detected load state
    Zl  = mean(Z_calc(m, ph));
    Zlt = mean(Z_true(m, ph));
    if isnan(real(r.Zth_true_final)), txt_true = '      n/a'; else, txt_true = sprintf('%8.4f %+8.4fj', real(r.Zth_true_final), imag(r.Zth_true_final)); end
    fprintf('  %d   | %8.4f %+8.4fj / %8.4f %+8.4fj          | %8.4f %+8.4fj / %s\n', ph, real(Zl), imag(Zl), real(Zlt), imag(Zlt), ...
            real(r.Zth_final), imag(r.Zth_final), txt_true);
end
fprintf('===========================================================================\n');

% --- Accuracy report: estimated vs actual Zth -----------------------------------------------------------------
f_zth = 48 + (52 - 48) * (res{1}.t_zth / T);        % frequency ramp of the test on the time grid of the graphs
fprintf('\n===========================================================================\n');
fprintf('  ACCURACY REPORT: ESTIMATED vs ACTUAL Zth\n');
fprintf('===========================================================================\n');
fprintf('Accuracy       = 100 * (1 - |final estimate - actual| / |actual|)\n');
fprintf('Max difference = largest |estimate - actual| over the whole graph, in Ohm. The reactance is compared at the same\n');
fprintf('                 frequency (X = 2*pi*f*L), so the ramp of the actual reactance does not count against the estimate.\n');
fprintf('Phase | Quantity |        Actual         |    Final estimate     | Accuracy | Max difference\n');
acc_summary = nan(n_ph, 6);                         % per phase: accuracy R, L, Zth [%] and max difference R, X, Zth [Ohm]
for ph = 1:n_ph
    r = res{ph};
    if isnan(Rth(ph)) || isnan(Lth(ph))
        fprintf('  %d   | no actual values given: no comparison\n', ph);
        continue;
    end
    d_t   = (real(r.Zth_calc) + 1i * 2*pi * f_zth .* r.L_trace) - (Rth(ph) + 1i * 2*pi * f_zth * Lth(ph));
    md_R  = max(abs(real(d_t)));                    % max difference of R_th [Ohm]
    md_X  = max(abs(imag(d_t)));                    % max difference of the reactance [Ohm]
    md_Z  = max(abs(d_t));                          % max difference of the complex Zth [Ohm]
    R_fin = real(r.Zth_final);
    L_fin = r.L_final;
    Z_act = Rth(ph) + 1i * 2*pi * f_nom * Lth(ph);  % Zth at the nominal frequency
    Z_fin = R_fin   + 1i * 2*pi * f_nom * L_fin;
    acc_R = max(0, 100 * (1 - abs(R_fin - Rth(ph)) / Rth(ph)));
    acc_L = max(0, 100 * (1 - abs(L_fin - Lth(ph)) / Lth(ph)));
    acc_Z = max(0, 100 * (1 - abs(Z_fin - Z_act) / abs(Z_act)));
    fprintf('  %d   | R_th     | %12.4f Ohm       | %12.4f Ohm       | %6.1f %% | %.4f Ohm (%.1f %% of actual)\n', ...
            ph, Rth(ph), R_fin, acc_R, md_R, 100 * md_R / Rth(ph));
    fprintf('      | L_th     | %14.4e H        | %14.4e H        | %6.1f %% | %.4f Ohm (%.1f %% of X at %g Hz)\n', ...
            Lth(ph), L_fin, acc_L, md_X, 100 * md_X / (2*pi*f_nom*Lth(ph)), f_nom);
    fprintf('      | Zth      | %8.4f %+8.4fj Ohm | %8.4f %+8.4fj Ohm | %6.1f %% | %.4f Ohm (%.1f %% of |Zth|)\n', ...
            real(Z_act), imag(Z_act), real(Z_fin), imag(Z_fin), acc_Z, md_Z, 100 * md_Z / abs(Z_act));
    acc_summary(ph, :) = [acc_R, acc_L, acc_Z, md_R, md_X, md_Z];
end
fprintf('Over the three phases | lowest accuracy: R_th %.1f %%, L_th %.1f %%, Zth %.1f %%\n', ...
        min(acc_summary(:, 1)), min(acc_summary(:, 2)), min(acc_summary(:, 3)));
fprintf('                      | largest max difference: R %.4f Ohm, X %.4f Ohm, Zth %.4f Ohm\n', ...
        max(acc_summary(:, 4)), max(acc_summary(:, 5)), max(acc_summary(:, 6)));
fprintf('===========================================================================\n');

% --- Continuous Zth graphs: start at the start value, grey bands mark the switch transients (line interpolated there)
figure('Name', 'Thevenin Impedance (continuous estimate), 3 phases', 'NumberTitle', 'off', ...
       'Position', [90 50 1150 820]);
part_name  = {'R_{th}', 'X_{th}'};
part_label = {'Resistance (\Omega)', 'Reactance (\Omega)'};
part_col   = {'b', 'r'};
for ph = 1:n_ph
    r          = res{ph};
    part_trace = {real(r.Zth_calc), imag(r.Zth_calc)};
    part_start = {real(Zth_start(ph)), imag(Zth_start(ph))};
    part_ev    = {real(r.Zth_event), imag(r.Zth_event)};
    part_use   = {r.use_R, r.use_X};
    part_true  = {real(Zth_true(:, ph)), imag(Zth_true(:, ph))};
    ok_ev      = ~isnan(r.Zth_event);

    for p = 1:2
        subplot(n_ph, 2, 2*(ph-1) + p);
        plot(r.t_zth, part_trace{p}, [part_col{p} '-'], 'LineWidth', 2.5, ...
             'DisplayName', ['Estimated ' part_name{p} sprintf(' (final %.4f \\Omega)', part_trace{p}(end))]);
        hold on;
        plot(0, part_start{p}, [part_col{p} 's'], 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'Clipping', 'off', ...
             'DisplayName', 'start value');
        shown = ok_ev & part_use{p};
        if any(shown)
            plot(r.t_event(shown), part_ev{p}(shown), [part_col{p} 'o'], 'MarkerFaceColor', part_col{p}, ...
                 'DisplayName', 'switch estimate (used)');
        end
        shown = ok_ev & ~part_use{p};
        if any(shown)
            plot(r.t_event(shown), part_ev{p}(shown), [part_col{p} 'x'], 'MarkerSize', 9, 'LineWidth', 1.5, ...
                 'DisplayName', 'switch estimate (not used)');
        end
        if ~isnan(Rth(ph)) && ~isnan(Lth(ph))
            plot(t_update, part_true{p}, 'k--', 'LineWidth', 2.0, 'DisplayName', ['Actual ' part_name{p}]);
        end

        % y-range: every plotted value, but never zoomed in tighter than half the start value
        vals = [part_trace{p}; part_ev{p}(ok_ev)];
        if ~isnan(Rth(ph)) && ~isnan(Lth(ph)), vals = [vals; part_true{p}]; end   % keep the actual line in view
        span = max(max(vals) - min(vals), 0.5 * abs(part_start{p}));
        ylim(0.5 * (max(vals) + min(vals)) + 0.65 * span * [-1 1]);

        % grey band over the unclean gap around every switch
        yl = ylim;
        for e = 1:r.n_ev
            hp = patch(r.t_gap(e, [1 2 2 1]), yl([1 1 2 2]), [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            uistack(hp, 'bottom');
            if e == 1
                set(hp, 'DisplayName', 'switch transient (interpolated)');
            else
                set(hp, 'HandleVisibility', 'off');
            end
        end
        hold off;
        grid on;
        ylabel(part_label{p});
        title(sprintf('Phase %d: Thevenin %s', ph, part_name{p}));
        if ph == 1, legend('Location', 'best', 'FontSize', 7); end
        if ph == n_ph, xlabel('Time (s)'); end
        xlim([0, T]);
    end
end

% --- Thevenin inductance in henries: L_th = X_th / (2*pi*f) at the frequency of each measurement ----------------
figure('Name', 'Thevenin Inductance L_th (H), 3 phases', 'NumberTitle', 'off', 'Position', [110 60 900 820]);
for ph = 1:n_ph
    r       = res{ph};
    L_start = imag(Zth_start(ph)) / (2*pi*f_nom);
    ok_ev   = ~isnan(r.L_event);
    subplot(n_ph, 1, ph);
    plot(r.t_zth, r.L_trace, 'm-', 'LineWidth', 2.5, ...
         'DisplayName', ['Estimated L_{th}' sprintf(' (final %.4g H)', r.L_trace(end))]);
    hold on;
    plot(0, L_start, 'ms', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'Clipping', 'off', 'DisplayName', 'start value');
    shown = ok_ev & r.use_X;
    if any(shown)
        plot(r.t_event(shown), r.L_event(shown), 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'switch estimate (used)');
    end
    shown = ok_ev & ~r.use_X;
    if any(shown)
        plot(r.t_event(shown), r.L_event(shown), 'mx', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'switch estimate (not used)');
    end
    vals = [r.L_trace; r.L_event(ok_ev)];
    if ~isnan(Lth(ph))
        plot([0 T], [Lth(ph) Lth(ph)], 'k--', 'LineWidth', 2.0, 'DisplayName', 'Actual L_{th}');
        vals = [vals; Lth(ph)];                      % keep the actual line in view
    end
    span = max(max(vals) - min(vals), 0.5 * abs(L_start));
    ylim(0.5 * (max(vals) + min(vals)) + 0.65 * span * [-1 1]);

    yl = ylim;                                       % grey band over the unclean gap around every switch
    for e = 1:r.n_ev
        hp = patch(r.t_gap(e, [1 2 2 1]), yl([1 1 2 2]), [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        uistack(hp, 'bottom');
        if e == 1
            set(hp, 'DisplayName', 'switch transient (interpolated)');
        else
            set(hp, 'HandleVisibility', 'off');
        end
    end
    hold off;
    grid on;
    ylabel('Inductance (H)');
    title(sprintf('Phase %d: Thevenin inductance L_{th}', ph));
    if ph == 1, legend('Location', 'best', 'FontSize', 7); end
    if ph == n_ph, xlabel('Time (s)'); end
    xlim([0, T]);
end

% --- Detection view: per-update change of R and X against the thresholds (for tuning)
figure('Name', 'Switch detection from Z_load, 3 phases', 'NumberTitle', 'off', 'Position', [120 40 1150 820]);
for ph = 1:n_ph
    r = res{ph};

    subplot(n_ph, 2, 2*ph - 1);
    plot(t_update(1:end-1), r.dR_step, 'b-', 'DisplayName', '\Delta R per update');
    hold on;
    plot([0 T], [ dR_thresh(ph)  dR_thresh(ph)], 'k--', 'DisplayName', 'threshold');
    plot([0 T], [-dR_thresh(ph) -dR_thresh(ph)], 'k--', 'HandleVisibility', 'off');
    plot(r.t_event, zeros(size(r.t_event)), 'rv', 'MarkerFaceColor', 'r', 'DisplayName', 'detected switch');
    hold off;
    grid on;
    ylabel('\Delta R (\Omega)');
    title(sprintf('Phase %d: change of R between consecutive windows', ph));
    if ph == 1, legend('Location', 'best', 'FontSize', 7); end
    if ph == n_ph, xlabel('Time (s)'); end
    xlim([0, T]);

    subplot(n_ph, 2, 2*ph);
    plot(t_update(1:end-1), r.dX_step, 'r-', 'DisplayName', '\Delta X per update');
    hold on;
    plot([0 T], [ dX_thresh(ph)  dX_thresh(ph)], 'k--', 'DisplayName', 'threshold');
    plot([0 T], [-dX_thresh(ph) -dX_thresh(ph)], 'k--', 'HandleVisibility', 'off');
    plot(r.t_event, zeros(size(r.t_event)), 'rv', 'MarkerFaceColor', 'r', 'DisplayName', 'detected switch');
    hold off;
    grid on;
    ylabel('\Delta X (\Omega)');
    title(sprintf('Phase %d: change of X between consecutive windows', ph));
    if ph == 1, legend('Location', 'best', 'FontSize', 7); end
    if ph == n_ph, xlabel('Time (s)'); end
    xlim([0, T]);
end

% =========================================================================
%  LOCAL HELPER FUNCTION: IMPEDANCE OF A SERIES R-L-C BRANCH (L = 0: no inductor, C = Inf: no capacitor)
% =========================================================================
function Zb = branch_z(R, L, C, w)
    Zb = R + 1i * w * L;
    if isfinite(C), Zb = Zb - 1i / (w * C); end
end

% =========================================================================
%  LOCAL HELPER FUNCTION: SWITCH DETECTION + THEVENIN IMPEDANCE OF ONE PHASE
% =========================================================================
function r = thevenin_from_load_changes(Z_calc, Ik_vec, Nw_vec, tc_vec, f_est, Tw_vec, t_update, Zth_true, Zth_start, cfg)
    % Inputs: Z_calc and Ik_vec of ONE phase (num_updates x 1) and the window quantities shared by the phases.
    % Zth_true (NaN when unknown) and Zth_start belong to the same phase. cfg holds the settings.
    num_updates = numel(t_update);
    T = cfg.T;

    % --- 1. Detect the switching events from the changes in Z_load ----------------
    % Real and imaginary part are tested separately (not the magnitude of the complex change);
    % an increase or a decrease of either part counts.
    dR_step = diff(real(Z_calc));
    dX_step = diff(imag(Z_calc));
    s_flag  = find(abs(dR_step) > cfg.dR_thresh | abs(dX_step) > cfg.dX_thresh);   % step s = window s -> window s+1

    n_merge = ceil(cfg.max_Tw / cfg.T_step);          % flagged steps within one window length belong to one switch
    if isempty(s_flag)
        ev_first = zeros(0, 1);
        ev_last  = zeros(0, 1);
    else
        gap_big  = diff(s_flag(:)) > n_merge;         % true where the next flagged step belongs to a new switch
        ev_first = s_flag([true; gap_big]);           % first window touched by each switch
        ev_last  = s_flag([gap_big; true]) + 1;       % last window touched by each switch
    end
    n_ev    = numel(ev_first);
    t_event = 0.5 * (t_update(ev_first) + Tw_vec(ev_first) + t_update(ev_last));   % estimated switching instants (s)

    % --- 2. Windows untouched by a switch, and the load state of each window ------
    is_clean = true(num_updates, 1);
    t_gap    = zeros(n_ev, 2);                        % [start end] (s) of the unclean gap around each switch
    for e = 1:n_ev
        t_gap(e, :) = [t_update(max(ev_first(e) - cfg.n_guard, 1)), t_update(ev_last(e)) + cfg.t_settle];
        is_clean(t_update >= t_gap(e, 1) & t_update < t_gap(e, 2)) = false;
    end
    state = 1 + sum((1:num_updates)' > ev_last', 2);  % 1 = before the first switch, 2 = after it, ...

    % --- 3. Current phasor from the coefficients, source rotation removed ---------
    % Ik_calc is referenced to t = 0 of the simulation and grows with the window length.
    % Both cancel in Z = V/I, but not in P = I1/I2 (two different windows), so refer every
    % window to its own centre and normalise by its length first.
    Iph = (Ik_vec ./ Nw_vec) .* exp(1i * 2 * pi * f_est .* tc_vec);

    % CHANGED: voltage phasor of every window, treated the same way (V_k = Z_k * I_k).
    % It is used ONLY to track the source phase; P is still formed from the currents.
    Vph = Iph .* Z_calc;

    % I_1 and I_2 must also be compared at the same source phase, whatever the grid frequency did in
    % between. The source phase is read from the phase of Vph in the clean stretches and carried across
    % the unclean gap at each switch by a smoothing spline. No frequency model is assumed.
    % CHANGED: the phase is tracked from V, not I. The angle of I relative to the source, -arg(Z_L + Zth),
    % changes with frequency, at a different rate in each load state, so the spline took that drift for
    % source movement and arg(P) came out a few mrad off at the switches. V = E*Z_L/(Z_L + Zth) stays
    % almost in phase with the source because Zth << Z_L, so its angle changes far less with frequency.
    %
    % 3a. Phase of Vph, unwrapped inside every clean stretch. f_est only supplies the whole-cycle count
    %     between consecutive windows: it sits on a 0.1 Hz grid, so integrating it directly drifts by
    %     several hundredths of a radian, far more than the phase of Vph itself is off. Each stretch is
    %     anchored to a linear reference phase (ph_ref) so that all the numbers stay small.
    ph_ref = 2 * pi * mean(f_est) * (tc_vec - tc_vec(1));
    U = nan(num_updates, 1);
    for s = 1:n_ev+1
        kk = find(is_clean & state == s);
        if numel(kk) < 2, continue; end
        dtc  = diff(tc_vec(kk));
        dphi = angle(Vph(kk(2:end)) ./ Vph(kk(1:end-1)));                       % CHANGED: Vph instead of Iph
        ncyc = round(0.5 * (f_est(kk(2:end)) + f_est(kk(1:end-1))) .* dtc - dphi / (2*pi));
        U(kk) = ph_ref(kk(1)) + angle(Vph(kk(1)) * exp(-1i * ph_ref(kk(1)))) + [0; cumsum(dphi + 2*pi*ncyc)];   % CHANGED
    end
    ok   = ~isnan(U);
    rows = find(ok);
    nr   = numel(rows);

    % 3b. Smoothing spline with one free offset per load state:
    %         minimise  sum_clean ( theta_k + c_state(k) - U_k )^2  +  lambda * sum ( 3rd divided difference of theta )^2
    %     The offsets c take up the phase step that every load change causes in the voltage, so theta never
    %     has to follow a load step: it follows the source only, and it is one continuous curve through every
    %     switching boundary. The third-derivative penalty leaves constant frequency and linear ramps
    %     unpenalised, and random or non-linear wander slower than f_smooth is followed. Across each unclean
    %     gap theta is interpolated with minimum third-derivative energy.
    %     lambda = 36*(2*pi*f_smooth)^-6 puts the half-power point of the smoother at f_smooth (36 = (3!)^2).
    theta_src    = ph_ref;                % fallback when there is nothing to fit
    fit_rms_mrad = NaN;
    if nr >= 4
        n      = num_updates;
        lambda = 36 * (2*pi*cfg.f_smooth)^(-6);
        Dd = speye(n);                    % 3rd divided difference on the real (slightly uneven) window-centre times
        for j = 1:3
            Dd = spdiags(1 ./ (tc_vec(1+j:end) - tc_vec(1:end-j)), 0, n-j, n-j) * (Dd(2:end, :) - Dd(1:end-1, :));
        end
        has   = accumarray(state(ok), 1, [n_ev+1, 1]) > 0;   % load states that have clean windows
        s_ref = find(has, 1);                                % the first of them is the phase reference (offset 0)
        s_off = find(has);
        s_off(s_off == s_ref) = [];                          % states with a free offset
        M   = [sparse((1:nr)', rows, 1, nr, n), sparse(double(state(rows) == s_off'))];   % [theta; offsets] -> U
        H   = M' * M;
        H(1:n, 1:n) = H(1:n, 1:n) + lambda * (Dd' * Dd);
        rhs = M' * (U(rows) - ph_ref(rows));
        z   = H \ rhs;
        z   = z + H \ (rhs - H * z);                         % one refinement step (the system is badly scaled)
        theta_src    = z(1:n) + ph_ref;                      % continuous source phase (up to a constant)
        fit_rms_mrad = 1e3 * sqrt(mean((M * z - (U(rows) - ph_ref(rows))).^2));
    end
    Idr = Iph .* exp(-1i * theta_src);                       % currents with the source rotation removed (P is still I1/I2)

    % --- 4. Apply the formula at every detected switch ---------------------------
    Zth_event      = complex(nan(n_ev, 1), nan(n_ev, 1));
    P_event        = complex(nan(n_ev, 1), nan(n_ev, 1));
    Zth_true_event = complex(nan(n_ev, 1), nan(n_ev, 1));
    sig_R_event    = nan(n_ev, 1);   % predicted error of Re(Zth) from each switch
    sig_X_event    = nan(n_ev, 1);   % predicted error of Im(Zth) from each switch
    f_event        = nan(n_ev, 1);   % frequency at which each switch measured Zth (Hz)

    for e = 1:n_ev
        k_pre  = find(is_clean & state == e,   cfg.n_avg, 'last');    % previous state, just before the switch
        k_post = find(is_clean & state == e+1, cfg.n_avg, 'first');   % current state, after the transient
        if numel(k_pre) < cfg.n_avg || numel(k_post) < cfg.n_avg, continue; end   % e.g. start-up: no clean state before it

        ZL1 = mean(Z_calc(k_post));                   % load after the switch
        ZL2 = mean(Z_calc(k_pre));                    % load before the switch
        if abs(real(ZL1 - ZL2)) <= cfg.dR_thresh && abs(imag(ZL1 - ZL2)) <= cfg.dX_thresh
            continue;                                 % no net load change: false alarm
        end

        P = mean(Idr(k_post)) / mean(Idr(k_pre));     % P = I_1 / I_2
        P_event(e)        = P;
        Zth_event(e)      = (ZL2 - P * ZL1) / (P - 1);
        Zth_true_event(e) = mean(Zth_true(k_post));
        f_event(e)        = mean(f_est(k_post));

        % error gain: dZth = err_gain * (amplitude error + j*phase error of P)
        err_gain = (ZL1 - ZL2) * P / (P - 1)^2;
        sig_R_event(e) = hypot(cfg.sig_amp * real(err_gain), cfg.sig_phase * imag(err_gain));
        sig_X_event(e) = hypot(cfg.sig_amp * imag(err_gain), cfg.sig_phase * real(err_gain));
    end

    % --- 4b. Update only the part of Zth that each switch measured reliably --------
    % Comparisons with NaN are false, so switches that were skipped above are never used.
    use_R = sig_R_event <= cfg.tol_rel * abs(Zth_event);
    use_X = sig_X_event <= cfg.tol_rel * abs(Zth_event);

    R_now = real(Zth_start);  X_now = imag(Zth_start);    % stay at the start value until a switch has measured that part accurately enough
    R_after = nan(n_ev, 1);   X_after = nan(n_ev, 1);      % estimate in force after each switch
    for e = 1:n_ev
        uR = use_R(1:e);                              % switches so far that measured R_th reliably
        uX = use_X(1:e);
        if any(uR)
            R_now = sum(real(Zth_event(uR)) ./ sig_R_event(uR).^2) / sum(1 ./ sig_R_event(uR).^2);
        end
        if any(uX)
            X_now = sum(imag(Zth_event(uX)) ./ sig_X_event(uX).^2) / sum(1 ./ sig_X_event(uX).^2);
        end
        R_after(e) = R_now;  X_after(e) = X_now;
    end
    Zth_final = complex(R_now, X_now);                % estimate from all detected switches
    if any(use_R), sig_R_final = 1 / sqrt(sum(1 ./ sig_R_event(use_R).^2)); else, sig_R_final = NaN; end
    if any(use_X), sig_X_final = 1 / sqrt(sum(1 ./ sig_X_event(use_X).^2)); else, sig_X_final = NaN; end
    Zth_true_final = mean(Zth_true(is_clean & state == n_ev+1));

    % Inductance view of the reactance: L = X / (2*pi*f) at the frequency of each measurement (L does not change with
    % the frequency ramp, X does). Same switches, same weights and same start value as the reactance X_th.
    L_event     = imag(Zth_event) ./ (2*pi*f_event);
    sig_L_event = sig_X_event ./ (2*pi*f_event);
    L_start     = imag(Zth_start) / (2*pi*cfg.f_nom);
    L_now       = L_start;
    L_after     = nan(n_ev, 1);
    for e = 1:n_ev
        uX = use_X(1:e);
        if any(uX)
            L_now = sum(L_event(uX) ./ sig_L_event(uX).^2) / sum(1 ./ sig_L_event(uX).^2);
        end
        L_after(e) = L_now;
    end
    L_final = L_now;
    if any(use_X), sig_L_final = 1 / sqrt(sum(1 ./ sig_L_event(use_X).^2)); else, sig_L_final = NaN; end

    % Continuous Zth trace on a fine time grid: the start value, flat while nothing new is known, and a
    % smoothstep hand-over to the new estimate across the unclean gap around each switch (no valid data there).
    t_zth   = (0 : 0.01 : T)';
    R_trace = real(Zth_start) * ones(size(t_zth));
    X_trace = imag(Zth_start) * ones(size(t_zth));
    L_trace = L_start * ones(size(t_zth));
    R_prev  = real(Zth_start);  X_prev = imag(Zth_start);  L_prev = L_start;
    for e = 1:n_ev
        w = min(max((t_zth - t_gap(e, 1)) / (t_gap(e, 2) - t_gap(e, 1)), 0), 1);
        s = w.^2 .* (3 - 2*w);                        % 0 before the gap, 1 after it, zero slope at both ends
        R_trace = R_trace + (R_after(e) - R_prev) * s;
        X_trace = X_trace + (X_after(e) - X_prev) * s;
        L_trace = L_trace + (L_after(e) - L_prev) * s;
        R_prev = R_after(e);  X_prev = X_after(e);  L_prev = L_after(e);
    end
    Zth_calc = complex(R_trace, X_trace);             % continuous Zth, evaluated on t_zth

    r = struct('n_ev', n_ev, 't_event', t_event, 't_gap', t_gap, 'is_clean', is_clean, 'state', state, ...
               'dR_step', dR_step, 'dX_step', dX_step, 'fit_rms_mrad', fit_rms_mrad, ...
               'Zth_event', Zth_event, 'P_event', P_event, 'Zth_true_event', Zth_true_event, ...
               'sig_R_event', sig_R_event, 'sig_X_event', sig_X_event, 'use_R', use_R, 'use_X', use_X, ...
               'Zth_final', Zth_final, 'sig_R_final', sig_R_final, 'sig_X_final', sig_X_final, ...
               'Zth_true_final', Zth_true_final, 't_zth', t_zth, 'Zth_calc', Zth_calc, ...
               'f_event', f_event, 'L_event', L_event, 'sig_L_event', sig_L_event, 'L_final', L_final, ...
               'sig_L_final', sig_L_final, 'L_trace', L_trace);
end

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
        n_idx  = (0:N_curr - 1)';

        A_next = zeros(N_curr, 1, 'like', sig);
        D_next = zeros(N_curr, 1, 'like', sig);

        for k = 0:7                                  % same sums, same order as the original loop over n, but for all n at once
            seg    = A_shifted(mod(2*n_idx + k, N_prev) + 1);
            A_next = A_next + LoD(k+1) * seg;
            D_next = D_next + HiD(k+1) * seg;
        end
        A_curr = A_next;
        D_coeff{j} = D_next;
    end
    A6 = A_curr;
    D6 = D_coeff{6};
    D5 = D_coeff{5};
end