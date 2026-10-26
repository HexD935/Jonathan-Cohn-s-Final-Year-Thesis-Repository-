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

% Preallocate storage for multi-coefficient calculations
f_est   = zeros(num_updates, 1);
Tw_vec  = zeros(num_updates, 1);
Z_calc  = zeros(num_updates, 1, 'like', 1i);

% NEW: storage for the Thevenin calculation (all coefficient-domain quantities)
Ik_vec  = zeros(num_updates, 1, 'like', 1i);   % I coefficient sum of every window
Nw_vec  = zeros(num_updates, 1);               % number of samples in every window
tc_vec  = zeros(num_updates, 1);               % exact centre time of every window (s)

% =========================================================================
%  CHANGED: DECOMPOSE V AND I ONCE OVER THE FULL RECORD (not gated per window)
%  Only the complex reference q is gated inside the loop. Every window's
%  phasor then matches the direct time-domain result, and the loop runs about
%  three times faster (one DWT per window instead of three).
% =========================================================================
[A6_Vf, D6_Vf, D5_Vf] = dwt_sym4_l6(Vload, LoD, HiD);
[A6_If, D6_If, D5_If] = dwt_sym4_l6(Iout,  LoD, HiD);

% =========================================================================
%  SLIDING WINDOW PROCESSING LOOP
% =========================================================================
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

    % 4. Construct window gate and apply it to the complex reference q only
    %    (CHANGED: V and I are no longer gated; they were decomposed whole before the loop)
    g = zeros(N, 1);
    g(idx_start:idx_end) = 1;

    q_k = g .* exp(-1i * 2 * pi * fk * t);

    % 5. Execute 6-Level DWT (sym4) of the gated reference to extract A6, D6, and D5
    [A6_q, D6_q, D5_q] = dwt_sym4_l6(q_k, LoD, HiD);

    % Multi-Coefficient Impedance Calculation (A6 + D6 + D5), using the full-record V and I coefficients
    Vk_calc = sum(A6_Vf .* A6_q) + sum(D6_Vf .* D6_q) + sum(D5_Vf .* D5_q);
    Ik_calc = sum(A6_If .* A6_q) + sum(D6_If .* D6_q) + sum(D5_If .* D5_q);
    Z_calc(k) = Vk_calc / Ik_calc;

    % NEW: keep what the Thevenin calculation needs
    Ik_vec(k) = Ik_calc;                             % current coefficient sum
    Nw_vec(k) = idx_end - idx_start + 1;             % samples in this window
    tc_vec(k) = (idx_start + idx_end - 2) / (2*Fs);  % centre time of this window (s)
end

R_calc = real(Z_calc);
X_calc = imag(Z_calc);

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
%  QUANTITATIVE ERROR METRICS
% =========================================================================
switch_times = [2.5, 5.0, 7.5];
is_steady = true(num_updates, 1);
for st = switch_times
    is_steady = is_steady & ~(t_update >= (st - 0.20) & t_update <= (st + 0.22));
end

calc_rmse = @(calc, true_val) sqrt(mean((calc - true_val).^2));

rmse_R_all = calc_rmse(R_calc, R_true);
rmse_X_all = calc_rmse(X_calc, X_true);
rmse_R_ss  = calc_rmse(R_calc(is_steady), R_true(is_steady));
rmse_X_ss  = calc_rmse(X_calc(is_steady), X_true(is_steady));

fprintf('\n===========================================================================\n');
fprintf('         MULTI-COEFFICIENT METHOD PERFORMANCE SUMMARY                       \n');
fprintf('===========================================================================\n');
fprintf('Resistance RMS Error  | Overall: %.6f Ohm | Steady-State: %.6f Ohm\n', rmse_R_all, rmse_R_ss);
fprintf('Reactance RMS Error   | Overall: %.6f Ohm | Steady-State: %.6f Ohm\n', rmse_X_all, rmse_X_ss);
fprintf('===========================================================================\n');

% =========================================================================
%  PLOTTING: ACTUAL VS CALCULATED (MULTI-COEFFICIENT METHOD)
% =========================================================================
figure('Name', 'Actual vs Calculated Impedance (Multi-Coefficient)', 'NumberTitle', 'off');

% Resistance Plot
subplot(2, 1, 1);
plot(t_update, R_true, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Actual R');
hold on;
plot(t_update, R_calc, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Calculated R (Multi-Coeff)');
hold off;
grid on;
ylabel('Resistance (\Omega)');
title('Resistance: Actual vs. Calculated');
legend('Location', 'best');
xlim([0, T]);

% Reactance Plot
subplot(2, 1, 2);
plot(t_update, X_true, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Actual X');
hold on;
plot(t_update, X_calc, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Calculated X (Multi-Coeff)');
hold off;
grid on;
xlabel('Time (s)');
ylabel('Reactance (\Omega)');
title('Reactance: Actual vs. Calculated');
legend('Location', 'best');
xlim([0, T]);

% =========================================================================
%  NEW: SWITCH DETECTION + THEVENIN IMPEDANCE (DWT-coefficient domain)
%
%      P   = I_1 / I_2                    1 = state after the transient
%      Zth = (Z_L2 - P*Z_L1) / (P - 1)    2 = state before the transient
%
%  The switching events are detected from Z_calc itself (switch_times is NOT
%  used here). Only coefficient-domain results of the sliding window
%  (Ik_calc, Z_calc) enter the formula - no waveform is rebuilt in the time domain.
% =========================================================================
% --- Settings --------------------------------------------------------------
% A switch is spread over one window length (~10 updates), so the change per update is only
% a fraction (~1/10) of the full load step. Set the thresholds above the steady-state
% update-to-update noise of Z_calc and below ~1/10 of the smallest load step to catch.
dR_thresh = 0.008;   % [Ohm per update] change of the REAL part between consecutive windows that flags a switch
dX_thresh = 0.008;   % [Ohm per update] change of the IMAGINARY part that flags a switch
n_guard   = 2;       % windows before the first flagged one that are discarded as well
t_settle  = 0.22;    % [s] windows start counting as "post-transient" this long after the switch has passed through
n_avg     = 2;       % windows averaged on each side of a switch (5 x 20 ms = 0.1 s)
% Source-phase smoother (Section 3): frequency wander slower than f_smooth is followed. Higher follows faster
% wander but lets more measurement noise into the phase that is bridged across each switch; lower is smoother.
f_smooth  = 1.5;     % [Hz]

% Which part of Zth can a switch measure? An error in P reaches Zth through err_gain = dZ*P/(P-1)^2 (dZ = ZL1 - ZL2).
% The dominant error of P is a few mrad in its PHASE (bridging the source phase across the switch); that lands in
% the part of Zth the load change moved (R for a resistive change, X for a reactive one, with loads near 45 deg),
% so such a switch measures the OTHER part well. Each switch therefore updates Re(Zth) or Im(Zth) only if the
% predicted error of that part is small; parts measured by several switches are averaged with weight 1/error^2.
sig_phase = 3e-3;    % [rad] expected error of arg(P)
sig_amp   = 1e-3;    % [-]   expected relative error of |P|
tol_rel   = 0.10;    % a switch may update R_th (X_th) only if its predicted error is below this fraction of |Zth|
Zth_start = 0.1 + 0.2i;   % [Ohm] assumed Zth before any switch has measured it: where the Zth graph starts

% Source (Thevenin) impedance used in the Simulink model: Zth = Rth + j*2*pi*f*Lth.
% Fill these in to get error figures; leave them NaN to skip the comparison.
Rth = 0.1;           % [Ohm]
Lth = 0.0006366;           % [H]
Zth_true = Rth + 1i * 2 * pi * f_true * Lth;

% --- 1. Detect the switching events from the changes in Z_load --------------
% Real and imaginary part are tested separately (not the magnitude of the complex change);
% an increase or a decrease of either part counts.
dR_step = diff(real(Z_calc));
dX_step = diff(imag(Z_calc));
s_flag  = find(abs(dR_step) > dR_thresh | abs(dX_step) > dX_thresh);   % step s = window s -> window s+1

n_merge = ceil(max_Tw / T_step);                  % flagged steps within one window length belong to one switch
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

% --- 2. Windows untouched by a switch, and the load state of each window -----
is_clean = true(num_updates, 1);
t_gap    = zeros(n_ev, 2);        % [start end] (s) of the unclean gap around each switch
for e = 1:n_ev
    t_gap(e, :) = [t_update(max(ev_first(e) - n_guard, 1)), t_update(ev_last(e)) + t_settle];
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
% changes with frequency (about 10 mrad per Hz in the first load state), and at a different rate in
% each load state, so the spline took that drift for source movement and arg(P) came out 2-3 mrad off
% at every switch. V = E*Z_L/(Z_L + Zth) stays almost in phase with the source because Zth << Z_L
% (about 0.04 mrad per Hz in the first load state).
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
    lambda = 36 * (2*pi*f_smooth)^(-6);
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

for e = 1:n_ev
    k_pre  = find(is_clean & state == e,   n_avg, 'last');    % previous state, just before the switch
    k_post = find(is_clean & state == e+1, n_avg, 'first');   % current state, after the transient
    if numel(k_pre) < n_avg || numel(k_post) < n_avg, continue; end   % e.g. start-up: no clean state before it

    ZL1 = mean(Z_calc(k_post));                       % load after the switch
    ZL2 = mean(Z_calc(k_pre));                        % load before the switch
    if abs(real(ZL1 - ZL2)) <= dR_thresh && abs(imag(ZL1 - ZL2)) <= dX_thresh
        continue;                                     % no net load change: false alarm
    end

    P = mean(Idr(k_post)) / mean(Idr(k_pre));         % P = I_1 / I_2
    P_event(e)        = P;
    Zth_event(e)      = (ZL2 - P * ZL1) / (P - 1);
    Zth_true_event(e) = mean(Zth_true(k_post));

    % error gain: dZth = err_gain * (amplitude error + j*phase error of P)
    err_gain = (ZL1 - ZL2) * P / (P - 1)^2;
    sig_R_event(e) = hypot(sig_amp * real(err_gain), sig_phase * imag(err_gain));
    sig_X_event(e) = hypot(sig_amp * imag(err_gain), sig_phase * real(err_gain));
end

% --- 4b. Update only the part of Zth that each switch measured reliably --------
% Comparisons with NaN are false, so switches that were skipped above are never used.
use_R = sig_R_event <= tol_rel * abs(Zth_event);
use_X = sig_X_event <= tol_rel * abs(Zth_event);

R_now = real(Zth_start);  X_now = imag(Zth_start);    % stay at the assumed start value until a switch has measured that part
R_after = nan(n_ev, 1);   X_after = nan(n_ev, 1);      % estimate in force after each switch
for e = 1:n_ev
    uR = use_R(1:e);                                  % switches so far that measured R_th reliably
    uX = use_X(1:e);
    if any(uR)
        R_now = sum(real(Zth_event(uR)) ./ sig_R_event(uR).^2) / sum(1 ./ sig_R_event(uR).^2);
    end
    if any(uX)
        X_now = sum(imag(Zth_event(uX)) ./ sig_X_event(uX).^2) / sum(1 ./ sig_X_event(uX).^2);
    end
    R_after(e) = R_now;  X_after(e) = X_now;
end
Zth_final = complex(R_now, X_now);                    % estimate from all detected switches

% Continuous Zth trace on a fine time grid: the assumed start value, flat while nothing new is known, and a
% smoothstep hand-over to the new estimate across the unclean gap around each switch (no valid data there).
t_zth   = (0 : 0.01 : T)';
R_trace = real(Zth_start) * ones(size(t_zth));
X_trace = imag(Zth_start) * ones(size(t_zth));
R_prev  = real(Zth_start);  X_prev = imag(Zth_start);
for e = 1:n_ev
    w = min(max((t_zth - t_gap(e, 1)) / (t_gap(e, 2) - t_gap(e, 1)), 0), 1);
    s = w.^2 .* (3 - 2*w);                            % 0 before the gap, 1 after it, zero slope at both ends
    R_trace = R_trace + (R_after(e) - R_prev) * s;
    X_trace = X_trace + (X_after(e) - X_prev) * s;
    R_prev = R_after(e);  X_prev = X_after(e);
end
Zth_calc = complex(R_trace, X_trace);                 % continuous Zth, evaluated on t_zth
if any(use_R), sig_R_final = 1 / sqrt(sum(1 ./ sig_R_event(use_R).^2)); else, sig_R_final = NaN; end
if any(use_X), sig_X_final = 1 / sqrt(sum(1 ./ sig_X_event(use_X).^2)); else, sig_X_final = NaN; end

% --- 5. Report ----------------------------------------------------------------
fprintf('\n===========================================================================\n');
fprintf('     SWITCH DETECTION + THEVENIN IMPEDANCE (P = I1/I2, DWT domain)          \n');
fprintf('===========================================================================\n');
fprintf('%d event(s) detected: |dR| > %.4g Ohm or |dX| > %.4g Ohm between consecutive windows\n', ...
        n_ev, dR_thresh, dX_thresh);
fprintf('Source-phase spline (phase tracked from V): f_smooth = %.2g Hz, fit residual %.2f mrad rms', f_smooth, fit_rms_mrad);   % CHANGED
if fit_rms_mrad > 5, fprintf(' (above ~5 mrad: try a higher f_smooth)'); end
fprintf('\n');
fprintf('A switch updates R_th / X_th only if its predicted error is below %.0f%% of |Zth|\n', 100 * tol_rel);
fprintf('The Zth graph starts from the assumed Zth = %.4f %+.4fj Ohm\n', real(Zth_start), imag(Zth_start));
used_txt = {'not used', 'used'};
for e = 1:n_ev
    [dt_min, j] = min(abs(switch_times - t_event(e)));   % known times only used to show the detection error
    if dt_min < 0.1
        fprintf('Event %d at %.3f s (actual %.1f s)', e, t_event(e), switch_times(j));
    else
        fprintf('Event %d at %.3f s (no known switch nearby)', e, t_event(e));
    end
    if isnan(Zth_event(e))
        fprintf(' | not used (no clean state before/after, or no net load change)\n');
        continue;
    end
    fprintf(' | |P-1| = %.3f\n', abs(P_event(e) - 1));
    fprintf('    R_th = %.4f +- %.4f (%s)    X_th = %.4f +- %.4f (%s)\n', ...
            real(Zth_event(e)), sig_R_event(e), used_txt{use_R(e) + 1}, ...
            imag(Zth_event(e)), sig_X_event(e), used_txt{use_X(e) + 1});
    if ~isnan(Rth) && ~isnan(Lth)
        fprintf('    true  = %.4f %+.4fj | actual error: R %+.4f, X %+.4f Ohm\n', ...
                real(Zth_true_event(e)), imag(Zth_true_event(e)), ...
                real(Zth_event(e) - Zth_true_event(e)), imag(Zth_event(e) - Zth_true_event(e)));
    end
end
if any(use_R), txt_R = sprintf('%d switch(es) +- %.4f', sum(use_R), sig_R_final); else, txt_R = 'assumed start value'; end
if any(use_X), txt_X = sprintf('%d switch(es) +- %.4f', sum(use_X), sig_X_final); else, txt_X = 'assumed start value'; end
fprintf('Combined Zth = %.4f %+.4fj Ohm   (R_th: %s | X_th: %s)\n', real(Zth_final), imag(Zth_final), txt_R, txt_X);
if ~isnan(Rth) && ~isnan(Lth)
    Zth_true_final = mean(Zth_true(is_clean & state == n_ev+1));
    fprintf('Actual   Zth = %.4f %+.4fj Ohm   (error: R %+.4f, X %+.4f)\n', real(Zth_true_final), imag(Zth_true_final), ...
            real(Zth_final - Zth_true_final), imag(Zth_final - Zth_true_final));
    % NEW: X_th = 2*pi*f*Lth rises with the frequency ramp, so the last load state is not a fair reference for a
    % value measured at earlier switches. Also compare with the actual Zth at the switches that were used,
    % weighted in the same way as the estimate.
    if any(use_R)
        R_ref = sum(real(Zth_true_event(use_R)) ./ sig_R_event(use_R).^2) / sum(1 ./ sig_R_event(use_R).^2);
    else
        R_ref = NaN;
    end
    if any(use_X)
        X_ref = sum(imag(Zth_true_event(use_X)) ./ sig_X_event(use_X).^2) / sum(1 ./ sig_X_event(use_X).^2);
    else
        X_ref = NaN;
    end
    fprintf('Actual Zth at the switches used = %.4f %+.4fj Ohm   (error: R %+.4f, X %+.4f)\n', ...
            R_ref, X_ref, real(Zth_final) - R_ref, imag(Zth_final) - X_ref);
end
fprintf('===========================================================================\n');

% Continuous Zth graph: starts at the assumed Zth, grey bands mark the switch transients (line interpolated there)
figure('Name', 'Thevenin Impedance (continuous estimate)', 'NumberTitle', 'off');
part_name  = {'R_{th}', 'X_{th}'};
part_label = {'Resistance (\Omega)', 'Reactance (\Omega)'};
part_col   = {'b', 'r'};
part_trace = {real(Zth_calc), imag(Zth_calc)};
part_start = {real(Zth_start), imag(Zth_start)};
part_ev    = {real(Zth_event), imag(Zth_event)};
part_use   = {use_R, use_X};
part_true  = {real(Zth_true), imag(Zth_true)};
ok_ev      = ~isnan(Zth_event);

for p = 1:2
    subplot(2, 1, p);
    plot(t_zth, part_trace{p}, [part_col{p} '-'], 'LineWidth', 2.5, ...
         'DisplayName', ['Estimated ' part_name{p} sprintf(' (final %.4f \\Omega)', part_trace{p}(end))]);
    hold on;
    plot(0, part_start{p}, [part_col{p} 's'], 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'Clipping', 'off', ...
         'DisplayName', 'assumed start');
    shown = ok_ev & part_use{p};
    if any(shown)
        plot(t_event(shown), part_ev{p}(shown), [part_col{p} 'o'], 'MarkerFaceColor', part_col{p}, ...
             'DisplayName', 'switch estimate (used)');
    end
    shown = ok_ev & ~part_use{p};
    if any(shown)
        plot(t_event(shown), part_ev{p}(shown), [part_col{p} 'x'], 'MarkerSize', 9, 'LineWidth', 1.5, ...
             'DisplayName', 'switch estimate (not used)');
    end
    if ~isnan(Rth) && ~isnan(Lth)
        plot(t_update, part_true{p}, 'k--', 'LineWidth', 2.0, 'DisplayName', ['Actual ' part_name{p}]);
    end

    % y-range: every plotted value, but never zoomed in tighter than half the assumed start value
    vals = [part_trace{p}; part_ev{p}(ok_ev)];
    span = max(max(vals) - min(vals), 0.5 * abs(part_start{p}));
    ylim(0.5 * (max(vals) + min(vals)) + 0.65 * span * [-1 1]);

    % grey band over the unclean gap around every switch
    yl = ylim;
    for e = 1:n_ev
        hp = patch(t_gap(e, [1 2 2 1]), yl([1 1 2 2]), [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
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
    legend('Location', 'best');
    xlim([0, T]);
end
subplot(2, 1, 1); title('Thevenin impedance: continuous estimate, starting from the assumed Zth');
subplot(2, 1, 2); xlabel('Time (s)');

% Detection view: per-update change of R and X against the thresholds (for tuning)
figure('Name', 'Switch detection from Z_load', 'NumberTitle', 'off');

subplot(2, 1, 1);
plot(t_update(1:end-1), dR_step, 'b-', 'DisplayName', '\Delta R per update');
hold on;
plot([0 T], [ dR_thresh  dR_thresh], 'k--', 'DisplayName', 'threshold');
plot([0 T], [-dR_thresh -dR_thresh], 'k--', 'HandleVisibility', 'off');
plot(t_event, zeros(size(t_event)), 'rv', 'MarkerFaceColor', 'r', 'DisplayName', 'detected switch');
hold off;
grid on;
ylabel('\Delta R (\Omega)');
title('Change of R between consecutive windows');
legend('Location', 'best');
xlim([0, T]);

subplot(2, 1, 2);
plot(t_update(1:end-1), dX_step, 'r-', 'DisplayName', '\Delta X per update');
hold on;
plot([0 T], [ dX_thresh  dX_thresh], 'k--', 'DisplayName', 'threshold');
plot([0 T], [-dX_thresh -dX_thresh], 'k--', 'HandleVisibility', 'off');
plot(t_event, zeros(size(t_event)), 'rv', 'MarkerFaceColor', 'r', 'DisplayName', 'detected switch');
hold off;
grid on;
xlabel('Time (s)');
ylabel('\Delta X (\Omega)');
title('Change of X between consecutive windows');
legend('Location', 'best');
xlim([0, T]);

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