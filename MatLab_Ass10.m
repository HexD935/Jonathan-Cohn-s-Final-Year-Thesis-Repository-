% =========================================================================
%  MULTI-COEFFICIENT (SYM4 DWT) IMPEDANCE ESTIMATION - THREE PHASES, NON-IDEAL NEUTRAL
%
%  Estimates, from the Simulink output out.simout.Data:
%    1. the load impedance Z_load of every phase (sliding-window DWT coefficient sums),
%    2. the Thevenin impedance Zth of every phase from its load switching events, with the
%       voltage drop in the neutral wire removed,
%    3. the impedance Zn of the neutral wire.
%  The three phases share one window and one frequency estimate per update; everything else is done phase by phase.
%
%  Input: out.simout.Data with the columns
%    V1 I1 V2 I2 V3 I3          Zn is estimated from the three phases (In = I1 + I2 + I3), or
%    V1 I1 V2 I2 V3 I3 Vn In    Zn = Vn/In is measured (Vn: voltage across the neutral impedance,
%                               In: current in the neutral wire).
%
%  Neutral wire. The three load currents return through one wire of impedance Zn, so for phase p
%        E_p = Zth_p*I_p + V_p + Vn ,     Vn = Zn*(I1 + I2 + I3)
%  (E_p: source EMF, V_p: voltage across the load).
%    - Z_load_p = V_p/I_p is not affected by the neutral; switch detection uses it.
%    - Zth is computed from Z_L = (V_p + Vn)/I_p, the impedance seen from the source terminals (load plus
%      neutral drop). The result is then the phase's own Thevenin impedance, whatever the other phases do.
%    - A switch in another phase steps this phase's current through Vn. Those switches are detected, the
%      windows they touch are left out, and the source phase is bridged across them with the EMF proxy
%      V_p + Vn + Zth_start*I_p, which does not step.
%    - Between the two states of a switch the grid frequency moves, and with it the reactance of Zth;
%      zth_two_states() therefore solves for the constants R_th and L_th.
%
%  Without neutral sensors, Zn is the value for which the source EMF V_p + Zth_p*I_p + Zn*In of every
%  phase is independent of the load. It is found by a Levenberg-Marquardt fit (fit_neutral_from_phases,
%  no toolbox needed) and requires load changes in at least two phases and constant source amplitudes;
%  otherwise a warning is printed and Zn = 0 is used (the Zth values then still contain the neutral path).
%
%  The Thevenin formula uses only the DWT-coefficient sums of the windows; no waveform is rebuilt in the
%  time domain.
% =========================================================================

%% INPUT DATA
Fs = 8192;                                % sampling frequency [Hz]
T  = 10;                                  % record length [s]
N  = Fs * T;                              % number of samples
t  = (0:N-1)' / Fs;                       % time vector [s]

n_ph  = 3;                                % number of phases
n_col = size(out.simout.Data, 2);
if n_col < 2 * n_ph
    error('out.simout.Data must have the columns V1 I1 V2 I2 V3 I3 (and, if logged, Vn In); found %d columns.', n_col);
end
if n_col == 2 * n_ph + 1
    error('out.simout.Data has 7 columns: log either the six columns V1 I1 V2 I2 V3 I3, or those and both neutral signals Vn In (8 columns).');
end
if n_col > 2 * n_ph + 2
    warning('out.simout.Data has %d columns: only the first eight are used (V1 I1 V2 I2 V3 I3 Vn In).', n_col);
end
neutral_sensors = (n_col >= 2 * n_ph + 2);   % true: Vn and In are logged (columns 7 and 8); false: only V1 I1 V2 I2 V3 I3
Vmeas = out.simout.Data(:, [1 3 5]);      % V1 V2 V3   (N x 3) as read by the voltage sensors
Iout  = out.simout.Data(:, [2 4 6]);      % I1 I2 I3   (N x 3)

% Neutral wire settings
volt_sensor_at = 'load';                  % 'load': voltage sensors across the loads (the Simulink model)
                                          % 'bus': load top to source neutral (the sensors also read the neutral drop)
rigid_sources  = 'auto';                  % Zn fit without neutral columns - assumption about the sources:
                                          % 'rigid': the three sources rotate together (one frequency waveform, fixed phase offsets)
                                          % 'amplitude': only the amplitude of each source is constant (no common rotation assumed)
                                          % 'auto': 'rigid', unless the fit residual shows that the sources do not rotate together
if ~(strcmpi(rigid_sources, 'auto') || strcmpi(rigid_sources, 'rigid') || strcmpi(rigid_sources, 'amplitude'))
    error('rigid_sources must be ''auto'', ''rigid'' or ''amplitude'' (it is ''%s'').', rigid_sources);
end
if neutral_sensors
    % Columns 7 and 8 are logged. The neutral current equals I1 + I2 + I3 (Kirchhoff), which identifies the current column
    % (whichever order the two were logged in) and fixes its sign. The sign of the voltage column is then fixed so that the
    % neutral impedance dissipates power (Re(Zn) > 0).
    Xn   = out.simout.Data(:, 7:8);
    s_ph = sum(Iout, 2);                  % I1 + I2 + I3
    g_n  = zeros(1, 2);                   % least-squares gain of each neutral column against I1 + I2 + I3
    nrm  = zeros(1, 2);                   % size of each column
    dev  = zeros(1, 2);                   % relative deviation from that fit (zero for the neutral current)
    for c = 1:2
        g_n(c) = (s_ph' * Xn(:, c)) / (s_ph' * s_ph);
        nrm(c) = norm(Xn(:, c));
        dev(c) = norm(Xn(:, c) - g_n(c) * s_ph) / max(nrm(c), eps);
    end
    dev_fit = dev;                        % kept for the message below
    dev(abs(abs(g_n) - 1) > 0.05 | nrm == 0) = Inf;   % the neutral current has gain +-1 against I1+I2+I3; a column that is zero is no current
    [dev_In, c_In] = min(dev);
    if dev_In > 0.05
        error(['Neither column 7 nor column 8 is the neutral current I1+I2+I3 (gain against it: %+.3f and %+.3f, should be +-1; ' ...
               'deviation from it: %.0f %% and %.0f %%): check the neutral current sensor.'], g_n(1), g_n(2), 100 * dev_fit(1), 100 * dev_fit(2));
    end
    In = Xn(:, c_In) * sign(g_n(c_In));   % neutral current, positive in the direction of I1 + I2 + I3
    Vn = Xn(:, 3 - c_In);                 % voltage across the neutral impedance ...
    sgn_V = sign(Vn' * In);
    if sgn_V == 0, sgn_V = 1; end
    Vn = Vn * sgn_V;                      % ... as a drop in the direction of the neutral current
    fprintf('Neutral columns: column %d = In (equals I1+I2+I3 to within %.2f %%, sign %+d), column %d = Vn (sign %+d)\n', ...
            6 + c_In, 100 * dev_In, sign(g_n(c_In)), 9 - c_In, sgn_V);
    if strcmpi(volt_sensor_at, 'bus')
        Vload = Vmeas - Vn;               % the sensors read the load voltage plus the neutral drop
    else
        Vload = Vmeas;                    % V1 V2 V3 across the loads
    end
else
    % No neutral columns: the coefficient sum of In = I1 + I2 + I3 is formed after the window loop (the sums are linear),
    % and Vn = Zn*In is reconstructed once Zn has been estimated (below).
    if strcmpi(volt_sensor_at, 'bus')
        error('volt_sensor_at = ''bus'' needs the logged columns Vn and In: without them the load voltage cannot be separated from the neutral drop.');
    end
    fprintf('No neutral columns logged (%d columns): In = I1+I2+I3, and Zn is estimated from the response of the three phases.\n', n_col);
    Vload = Vmeas;                        % V1 V2 V3 across the loads
end

%% SYM4 FILTERS
% Decomposition filters of the sym4 wavelet: low-pass (LoD) and high-pass (HiD)
LoD = [ ...
   -0.07576571478927333, -0.02963552764599851,  0.49761866763201545,  0.80373875180591610, ...
    0.29785779560527736, -0.09921954357684722, -0.01260396726203783,  0.03222310060404270 ];
HiD = [ ...
   -0.03222310060404270, -0.01260396726203783,  0.09921954357684722,  0.29785779560527736, ...
   -0.80373875180591610,  0.49761866763201545,  0.02963552764599851, -0.07576571478927333 ];

%% SLIDING-WINDOW SETUP
% Adaptive 10-cycle window (length 10/fk), updated every 20 ms
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
Vk_vec  = zeros(num_updates, n_ph, 'like', 1i);   % load-voltage coefficient sum of every phase
Vnk_vec = zeros(num_updates, 1, 'like', 1i);      % coefficient sum of the voltage across the neutral impedance
Ink_vec = zeros(num_updates, 1, 'like', 1i);      % coefficient sum of the neutral current
Zn_calc = zeros(num_updates, 1, 'like', 1i);      % neutral impedance Vn/In

% DWT of V and I of every phase (and of Vn, In when logged), computed once over the full record.
% Only the complex reference q is gated inside the window loop, so each window needs a single DWT.
A6_Vf = zeros(N/2^6, n_ph);  D6_Vf = zeros(N/2^6, n_ph);  D5_Vf = zeros(N/2^5, n_ph);
A6_If = zeros(N/2^6, n_ph);  D6_If = zeros(N/2^6, n_ph);  D5_If = zeros(N/2^5, n_ph);
for ph = 1:n_ph
    [A6_Vf(:, ph), D6_Vf(:, ph), D5_Vf(:, ph)] = dwt_sym4_l6(Vload(:, ph), LoD, HiD);
    [A6_If(:, ph), D6_If(:, ph), D5_If(:, ph)] = dwt_sym4_l6(Iout(:, ph),  LoD, HiD);
end
if neutral_sensors
    [A6_Vn_f, D6_Vn_f, D5_Vn_f] = dwt_sym4_l6(Vn, LoD, HiD);
    [A6_In_f, D6_In_f, D5_In_f] = dwt_sym4_l6(In, LoD, HiD);
end

%% SLIDING-WINDOW LOOP
for k = 1:num_updates
    t_start = t_update(k);
    idx_start = round(t_start * Fs) + 1;

    % 1. Fundamental frequency over a nominal 200 ms interval (summed spectral energy of the three phase voltages)
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

    % 4. Window gate and complex reference q (the same for the three phases)
    g = zeros(N, 1);
    g(idx_start:idx_end) = 1;
    Nw_k = idx_end - idx_start + 1;                  % samples in this window
    q_k = g .* exp(-1i * 2 * pi * fk * t);
    [A6_q, D6_q, D5_q] = dwt_sym4_l6(q_k, LoD, HiD);

    % 5. A6 + D6 + D5 coefficient sums of every phase
    for ph = 1:n_ph
        Vk_calc = sum(A6_Vf(:, ph) .* A6_q) + sum(D6_Vf(:, ph) .* D6_q) + sum(D5_Vf(:, ph) .* D5_q);
        Ik_calc = sum(A6_If(:, ph) .* A6_q) + sum(D6_If(:, ph) .* D6_q) + sum(D5_If(:, ph) .* D5_q);
        Z_calc(k, ph) = Vk_calc / Ik_calc;
        Vk_vec(k, ph) = Vk_calc;                     % kept for the Thevenin calculation
        Ik_vec(k, ph) = Ik_calc;
    end

    % Same sums for the voltage across the neutral impedance and the neutral current (when logged)
    if neutral_sensors
        Vnk_vec(k) = sum(A6_Vn_f .* A6_q) + sum(D6_Vn_f .* D6_q) + sum(D5_Vn_f .* D5_q);
        Ink_vec(k) = sum(A6_In_f .* A6_q) + sum(D6_In_f .* D6_q) + sum(D5_In_f .* D5_q);
        Zn_calc(k) = Vnk_vec(k) / Ink_vec(k);
    end

    Nw_vec(k) = Nw_k;                                % samples in this window
    tc_vec(k) = (idx_start + idx_end - 2) / (2*Fs);  % centre time of this window (s)
end

if ~neutral_sensors
    Ink_vec = sum(Ik_vec, 2);                        % coefficient sum of In = I1 + I2 + I3 (the sums are linear)
end

R_calc = real(Z_calc);                               % num_updates x n_ph
X_calc = imag(Z_calc);

%% GROUND TRUTH (comparison only; not used by the estimator)
t_mid  = t_update + (Tw_vec / 2);
f_true = 48 + (52 - 48) * (t_mid / T);

% Local solver of the Simulink model (Solver Configuration block, 'Use local solver'). A fixed-step solver replaces s = j*w
% by a discrete approximation, so the simulated circuit differs slightly from the ideal one: under Backward Euler every
% inductor gains a series resistance (2*pi*f)^2*h*L/2 and every capacitor one of h/(2*C); the trapezoidal rule adds none.
% The estimator measures the simulated circuit and needs no correction. Z_true is the impedance of the ideal loads and
% Z_solv that of the same loads under the solver (both for the comparison only).
solver_type = 'BE';       % 'BE': Backward Euler, 'TR': Trapezoidal rule, 'none': exact (the model uses a continuous solver)
solver_h    = 1/8192;     % [s] step of the local solver

% Loads (comparison only): four series R-L-C branches per phase, connected in parallel; branch 1 is always on, branches 2, 3
% and 4 are added at the times in switch_times (one row per phase). Columns = [phase 1, phase 2, phase 3];
% L = 0: no inductor; C = Inf: no capacitor.
switch_times = [2.5  5.0  7.5;                    % phase 1
                2.0  6.0  8.0;                    % phase 2
                2.0  3.0  7.0];                   % phase 3
%        phase 1    phase 2   phase 3
R1 = [   1          0.8       0.6      ];   L1 = [ 0.00318   0         0.00191  ];   C1 = [ Inf       Inf       Inf      ];
R2 = [   10         0.2       1.5      ];   L2 = [ 0         0.00318   0        ];   C2 = [ Inf       Inf       Inf      ];
R3 = [   2          0.5       0.1      ];   L3 = [ 0.0318    0         0        ];   C3 = [ Inf       0.006366  0.003979 ];
R4 = [   5          2.0       0.3      ];   L4 = [ 0         0         0.001273 ];   C4 = [ 0.01592   Inf       Inf      ];
Z_true = zeros(num_updates, n_ph, 'like', 1i);
Z_solv = zeros(num_updates, n_ph, 'like', 1i);    % the same loads as the local solver implements them

for ph = 1:n_ph
    for k = 1:num_updates
        w  = 2 * pi * f_true(k);
        Zb = [branch_z(R1(ph), L1(ph), C1(ph), w), ...
              branch_z(R2(ph), L2(ph), C2(ph), w), ...
              branch_z(R3(ph), L3(ph), C3(ph), w), ...
              branch_z(R4(ph), L4(ph), C4(ph), w)];
        Zs = [branch_z_solver(R1(ph), L1(ph), C1(ph), w, solver_h, solver_type), ...
              branch_z_solver(R2(ph), L2(ph), C2(ph), w, solver_h, solver_type), ...
              branch_z_solver(R3(ph), L3(ph), C3(ph), w, solver_h, solver_type), ...
              branch_z_solver(R4(ph), L4(ph), C4(ph), w, solver_h, solver_type)];
        n_on = 1 + sum(t_mid(k) >= switch_times(ph, :));   % branches connected at this time
        Z_true(k, ph) = 1 / sum(1 ./ Zb(1:n_on));    % parallel combination
        Z_solv(k, ph) = 1 / sum(1 ./ Zs(1:n_on));
    end
end

R_true = real(Z_true);
X_true = imag(Z_true);
R_solv = real(Z_solv);
X_solv = imag(Z_solv);

%% ERROR OF Z_LOAD (per phase)
is_steady = true(num_updates, n_ph);              % steady-state windows of every phase (its own switching times)
for ph = 1:n_ph
    for st = switch_times(ph, :)
        is_steady(:, ph) = is_steady(:, ph) & ~(t_update >= (st - 0.20) & t_update <= (st + 0.22));
    end
end

calc_rmse = @(calc, true_val) sqrt(mean((calc - true_val).^2));
rmse_R_all = zeros(1, n_ph);  rmse_X_all = zeros(1, n_ph);
rmse_R_ss  = zeros(1, n_ph);  rmse_X_ss  = zeros(1, n_ph);
rmse_R_solv = zeros(1, n_ph);  rmse_X_solv = zeros(1, n_ph);   % steady-state error against the solver-equivalent impedance
off_R_solv  = zeros(1, n_ph);  off_X_solv  = zeros(1, n_ph);   % mean shift (solver minus ideal) in the steady-state windows
for ph = 1:n_ph
    ss = is_steady(:, ph);
    rmse_R_all(ph) = calc_rmse(R_calc(:, ph), R_true(:, ph));
    rmse_X_all(ph) = calc_rmse(X_calc(:, ph), X_true(:, ph));
    rmse_R_ss(ph)  = calc_rmse(R_calc(ss, ph), R_true(ss, ph));
    rmse_X_ss(ph)  = calc_rmse(X_calc(ss, ph), X_true(ss, ph));
    rmse_R_solv(ph) = calc_rmse(R_calc(ss, ph), R_solv(ss, ph));
    rmse_X_solv(ph) = calc_rmse(X_calc(ss, ph), X_solv(ss, ph));
    off_R_solv(ph)  = mean(R_solv(ss, ph) - R_true(ss, ph));
    off_X_solv(ph)  = mean(X_solv(ss, ph) - X_true(ss, ph));
end

fprintf('\n===========================================================================\n');
fprintf('         MULTI-COEFFICIENT METHOD PERFORMANCE SUMMARY (3 PHASES)            \n');
fprintf('===========================================================================\n');
for ph = 1:n_ph
    fprintf('Phase %d | Resistance RMS Error | Overall: %.6f Ohm | Steady-State: %.6f Ohm\n', ph, rmse_R_all(ph), rmse_R_ss(ph));
    fprintf('        | Reactance RMS Error  | Overall: %.6f Ohm | Steady-State: %.6f Ohm\n', rmse_X_all(ph), rmse_X_ss(ph));
end
if ~strcmpi(solver_type, 'none')
    fprintf('---------------------------------------------------------------------------\n');
    fprintf('Local solver of the model: %s, step %g s. The errors above are against the IDEAL loads, but the solver\n', solver_type, solver_h);
    fprintf('changes the loads that are simulated, and the estimator measures the simulated ones. Steady-state error\n');
    fprintf('against the impedance the solver implements (what the method itself gets wrong), and the shift it causes:\n');
    for ph = 1:n_ph
        fprintf('Phase %d | R: %.6f Ohm | X: %.6f Ohm | shift (solver minus ideal, mean): R %+.6f Ohm, X %+.6f Ohm\n', ...
                ph, rmse_R_solv(ph), rmse_X_solv(ph), off_R_solv(ph), off_X_solv(ph));
    end
end
fprintf('===========================================================================\n');

%% FIGURE: Z_LOAD, ACTUAL VS CALCULATED
figure('Name', 'Actual vs Calculated Impedance (Multi-Coefficient), 3 phases', 'NumberTitle', 'off', ...
       'Position', [60 60 1150 820]);
for ph = 1:n_ph
    % Resistance
    subplot(n_ph, 2, 2*ph - 1);
    plot(t_update, R_true(:, ph), 'k--', 'LineWidth', 2.0, 'DisplayName', 'Actual R');
    hold on;
    if ~strcmpi(solver_type, 'none')
        plot(t_update, R_solv(:, ph), ':', 'Color', [0 0.6 0], 'LineWidth', 2.0, 'DisplayName', ['Actual R under ' solver_type ' solver']);
    end
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
    if ~strcmpi(solver_type, 'none')
        plot(t_update, X_solv(:, ph), ':', 'Color', [0 0.6 0], 'LineWidth', 2.0, 'DisplayName', ['Actual X under ' solver_type ' solver']);
    end
    plot(t_update, X_calc(:, ph), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Calculated X (Multi-Coeff)');
    hold off;
    grid on;
    ylabel('Reactance (\Omega)');
    title(sprintf('Phase %d: Reactance, actual vs. calculated', ph));
    if ph == 1, legend('Location', 'best', 'FontSize', 8); end
    if ph == n_ph, xlabel('Time (s)'); end
    xlim([0, T]);
end

%% SWITCH DETECTION AND THEVENIN IMPEDANCE OF EACH PHASE
%  P   = I1/I2                      (1: state after the transient, 2: state before it)
%  Zth = (Z_L2 - P*Z_L1)/(P - 1)
%  with Z_L = (V_load + Vn)/I, the load plus the neutral drop as seen from the source terminals (see the header).
%
%  The switches of every phase are detected from its own Z_calc by detect_switches() (switch_times is not used for
%  that). thevenin_from_load_changes() then processes each phase; it also receives the switches of the other two
%  phases, because they act on this phase through the neutral. Only coefficient-domain results (the Vk, Ik and Vnk
%  sums) enter the formula; no waveform is rebuilt in the time domain.
% Settings (the same for the three phases, except where marked per phase)
% Detection thresholds per phase [phase 1, phase 2, phase 3]: a switch is spread over one window length (~10 updates), so the
% change per update is only ~1/10 of the load step. Set each threshold above the steady-state update-to-update noise of
% Z_calc of that phase and below ~1/10 of its smallest load step.
dR_thresh = [0.008 0.004 0.004];   % [Ohm per update] change of the real part between consecutive windows that flags a switch
dX_thresh = [0.008 0.004 0.004];   % [Ohm per update] change of the imaginary part that flags a switch
n_guard   = 2;       % windows before the first flagged one that are discarded as well
t_settle  = 0.22;    % [s] windows start counting as "post-transient" this long after the switch has passed through
n_avg     = 5;       % windows averaged on each side of a switch (5 x 20 ms = 0.1 s)
n_min     = 3;       % fewest clean windows accepted on either side of a switch (up to n_avg are used): the switches of the
                     % other phases shorten the clean stretches between this phase's own switches
t_settle_f = 0.10;   % [s] windows are also left out this long after a switch in another phase. That switch reaches this phase
                     % only through the neutral (a step of a few percent in its current), so less waiting is needed than after
                     % this phase's own switches (t_settle)
% Source-phase smoother: frequency wander slower than f_smooth is followed. A higher value follows faster wander but lets
% more measurement noise into the phase bridged across each switch; a lower value is smoother. 0.75 Hz suits the linear
% ramp; for random frequency wander use 1.5 - 3 Hz as far as its speed requires. The printed fit residual shows when the
% smoother is too stiff (above ~5 mrad).
f_smooth  = 0.75;    % [Hz]

% Reliability gate. An error in P reaches Zth through err_gain = dZ*P/(P-1)^2 (dZ = ZL1 - ZL2). The dominant error of P is
% in its phase (bridging the source phase across the switch); it lands in the part of Zth that the load change moved (R for
% a resistive change, X for a reactive one, with loads near 45 deg), so such a switch measures the other part well. Each
% switch therefore updates Re(Zth) or Im(Zth) only if the predicted error of that part is small; parts measured by several
% switches are averaged with weight 1/error^2. The two sigmas below are conservative estimates of the error of P.
sig_phase = 3e-3;    % [rad] expected error of arg(P)
sig_amp   = 1e-3;    % [-]   expected relative error of |P|
tol_rel   = 0.05;    % accuracy required for a Zth change: a switch may update R_th (X_th) only if its predicted error
                     % is below this fraction of |Zth| (smaller = stricter)

% Actual source (Thevenin) impedance of every phase [phase 1, phase 2, phase 3], Zth = Rth + j*2*pi*f*Lth. Used for the
% dashed 'Actual' lines, the error figures and as the start value of every Zth graph (below).
Rth = [0.1, 0.03, 0.05];                            % [Ohm]
Lth = [0.2/(2*pi*50), 0.000127, 0.000191];          % [H]   (phase 1: 0.2 Ohm of reactance at 50 Hz)
f_nom     = 50;                                     % [Hz] nominal frequency: start value, and X <-> L conversion of the start value
Zth_start = Rth + 1i * 2 * pi * f_nom * Lth;        % [Ohm] where each Zth graph starts: the actual value at f_nom;
                                                    % a switch replaces it only when it measures that part accurately enough
Zth_true = Rth + 1i * 2 * pi * (f_true .* Lth);   % num_updates x n_ph

% Actual impedance of the neutral wire (comparison only, not used by the estimator): 0.05 Ohm + 0.05 Ohm of reactance at 50 Hz.
Rn_act = 0.05;                                      % [Ohm]
Ln_act = 0.05 / (2*pi*50);                          % [H]

% Settings passed to the helper functions
cfg = struct('dR_thresh', dR_thresh(1), 'dX_thresh', dX_thresh(1), 'n_guard', n_guard, 't_settle', t_settle, ...
             't_settle_f', t_settle_f, 'n_avg', n_avg, 'n_min', n_min, 'f_smooth', f_smooth, ...
             'sig_phase', sig_phase, 'sig_amp', sig_amp, 'tol_rel', tol_rel, 'f_nom', f_nom, ...
             'T', T, 'T_step', T_step, 'max_Tw', max_Tw);

% Switch detection for every phase, from its own Z_load = V/I (not affected by the neutral)
ev_first = cell(n_ph, 1);  ev_last = cell(n_ph, 1);    % first / last window touched by each switch
dR_step  = cell(n_ph, 1);  dX_step = cell(n_ph, 1);    % change of R and X between consecutive windows
for ph = 1:n_ph
    cfg.dR_thresh = dR_thresh(ph);                    % detection thresholds of this phase
    cfg.dX_thresh = dX_thresh(ph);
    [ev_first{ph}, ev_last{ph}, dR_step{ph}, dX_step{ph}] = detect_switches(Z_calc(:, ph), cfg);
end

% Without neutral columns: estimate Zn from the three phases and reconstruct Vn = Zn*In
if ~neutral_sensors
    rot_all = exp(1i * 2 * pi * f_est .* tc_vec);     % window phasors referred to the window centres (as in the Thevenin function)
    Vph_all = (Vk_vec ./ Nw_vec) .* rot_all;
    Iph_all = (Ik_vec ./ Nw_vec) .* rot_all;
    ok_fit  = all(isfinite(Vph_all) & isfinite(Iph_all), 2) & any(abs(Vph_all) > 0, 2);   % windows with usable phasors
    st_all  = ones(num_updates, n_ph);                % load state of every window, phase by phase (1 = before its first switch)
    for ph = 1:n_ph
        st_all(:, ph) = 1 + sum((1:num_updates)' > ev_last{ph}', 2);
        for e = 1:numel(ev_first{ph})
            g0 = t_update(max(ev_first{ph}(e) - n_guard, 1));
            g1 = t_update(ev_last{ph}(e)) + t_settle;
            ok_fit(t_update >= g0 & t_update < g1) = false;   % leave out the transients of the switches of all phases
        end
    end
    nfit    = fit_neutral_from_phases(Vph_all(ok_fit, :), Iph_all(ok_fit, :), f_est(ok_fit), st_all(ok_fit, :), rigid_sources);
    Zn_calc = nfit.Rn + 1i * nfit.Xn * f_est / 50;    % Zn at the frequency of each window (reactance proportional to f)
    Vnk_vec = Zn_calc .* Ink_vec;                     % reconstructed coefficient sum of Vn
end

% Thevenin impedance of every phase from its own switches and those of the other two phases
res = cell(n_ph, 1);
for ph = 1:n_ph
    other = setdiff(1:n_ph, ph);                      % the other two phases
    cfg.dR_thresh = dR_thresh(ph);
    cfg.dX_thresh = dX_thresh(ph);
    res{ph} = thevenin_from_load_changes(Vk_vec(:, ph), Vnk_vec, Ik_vec(:, ph), Nw_vec, tc_vec, f_est, Tw_vec, t_update, ...
                                         ev_first{ph}, ev_last{ph}, vertcat(ev_first{other}), vertcat(ev_last{other}), ...
                                         Zth_true(:, ph), Zth_start(ph), cfg);
end

% Neutral wire impedance Zn: value of every window (Zn_calc) and its estimate, the mean over the windows untouched by a switch
clean_all = true(num_updates, 1);                     % windows untouched by a switch in any phase
for ph = 1:n_ph
    clean_all = clean_all & res{ph}.is_clean;
end
if ~any(clean_all)
    warning('Every window is within a switch transient: the neutral impedance is averaged over all windows instead.');
    clean_all = true(num_updates, 1);
end
Rn_calc = real(Zn_calc);
Ln_calc = imag(Zn_calc) ./ (2 * pi * f_est);          % L = X/(2*pi*f) at the frequency of each window
Rn_est  = mean(Rn_calc(clean_all));                   % estimate: mean over the switch-free windows
Ln_est  = mean(Ln_calc(clean_all));
Zn_true = Rn_act + 1i * 2 * pi * (f_true * Ln_act);   % actual Zn at the frequency of each window (comparison only)
Zn_err  = Zn_calc - Zn_true;
Zn_nom  = Rn_act + 1i * 2*pi * f_nom * Ln_act;        % actual Zn at the nominal frequency
Zn_fin  = Rn_est + 1i * 2*pi * f_nom * Ln_est;        % estimated Zn at the nominal frequency
acc_Rn  = max(0, 100 * (1 - abs(Rn_est - Rn_act) / Rn_act));
acc_Ln  = max(0, 100 * (1 - abs(Ln_est - Ln_act) / Ln_act));
acc_Zn  = max(0, 100 * (1 - abs(Zn_fin - Zn_nom) / abs(Zn_nom)));

%% REPORTS
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
    fprintf('Source-phase spline (phase tracked from the EMF proxy V + Vn + Zth_start*I): f_smooth = %.2g Hz, fit residual %.2f mrad rms', f_smooth, r.fit_rms_mrad);
    if r.fit_rms_mrad > 5, fprintf(' (above ~5 mrad: try a higher f_smooth)'); end
    fprintf('\n');
    fprintf('The Zth graph starts from the start value %.4f %+.4fj Ohm\n', real(Zth_start(ph)), imag(Zth_start(ph)));
    if isempty(r.t_foreign)
        fprintf('No switches detected in the other phases.\n');
    else
        fprintf('Switches in the other phases (they reach this phase through the neutral; their windows are left out): t = %s s\n', ...
                sprintf('%.2f ', r.t_foreign));
    end
    for e = 1:r.n_ev
        [dt_min, j] = min(abs(switch_times(ph, :) - r.t_event(e)));   % known times only used to show the detection error
        if dt_min < 0.1
            fprintf('Event %d at %.3f s (actual %.1f s)', e, r.t_event(e), switch_times(ph, j));
        else
            fprintf('Event %d at %.3f s (no known switch nearby)', e, r.t_event(e));
        end
        if isnan(r.Zth_event(e))
            fprintf(' | not used (clean windows before/after: %d/%d, at least %d needed; or no net load change)\n', ...
                    r.n_pre(e), r.n_post(e), n_min);
            continue;
        end
        fprintf(' | |P-1| = %.3f\n', abs(r.P_event(e) - 1));
        fprintf('    R_th = %.4f +- %.4f [limit %.4f] (%s)    X_th = %.4f +- %.4f [limit %.4f] (%s)\n', ...
                real(r.Zth_event(e)), r.sig_R_event(e), tol_rel * abs(r.Zth_event(e)), used_txt{r.use_R(e) + 1}, ...
                imag(r.Zth_event(e)), r.sig_X_event(e), tol_rel * abs(r.Zth_event(e)), used_txt{r.use_X(e) + 1});
        fprintf('    L_th = %.4e H  (X_th at %.2f Hz%s)\n', r.L_event(e), r.f_event(e), ...
                sprintf('; actual %.4e H', Lth(ph)));
        fprintf('    same switch WITHOUT the neutral correction (Z_L = V/I): R_th = %.4f, X_th = %.4f  (neutral path included)\n', ...
                real(r.Zth_blind(e)), imag(r.Zth_blind(e)));
        fprintf('    plain formula, Zth taken equal in both states: R_th = %.4f, X_th = %.4f  (the grid moved %+.2f Hz between them)\n', ...
                real(r.Zth_plain(e)), imag(r.Zth_plain(e)), r.df_event(e));
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
fprintf('Neutral wire Zn at %g Hz [Ohm]: estimate %8.4f %+8.4fj / actual %8.4f %+8.4fj\n', f_nom, ...
        real(Zn_fin), imag(Zn_fin), real(Zn_nom), imag(Zn_nom));
fprintf('===========================================================================\n');

% Accuracy report: estimated vs actual Zth
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
if strcmpi(solver_type, 'BE')
    dR_solv_th = (2*pi*f_nom)^2 * solver_h * Lth / 2;      % series resistance the solver adds to each source inductance [Ohm]
    dR_solv_n  = (2*pi*f_nom)^2 * solver_h * Ln_act / 2;   % ... and to the neutral inductance [Ohm]
    fprintf('Solver note: under Backward Euler (step %g s) the simulated source inductances gain a series resistance\n', solver_h);
    fprintf('(2*pi*f)^2*h*L/2: at %g Hz %.4f / %.4f / %.4f Ohm in phases 1 / 2 / 3 (neutral wire %.4f Ohm). That part of\n', ...
            f_nom, dR_solv_th(1), dR_solv_th(2), dR_solv_th(3), dR_solv_n);
    fprintf('the R_th and R_n differences comes from the solver, not from the estimator. The reactances shrink by only %.3f %%.\n', ...
            100 * (2*pi*f_nom*solver_h)^2 / 6);
    fprintf('===========================================================================\n');
end

% Neutral wire impedance: measured as Vn/In, or estimated from the three phases when Vn and In are not logged
fprintf('\n===========================================================================\n');
if neutral_sensors
    fprintf('  NEUTRAL WIRE IMPEDANCE Zn = Vn/In  (%d of %d windows are free of any switch)\n', sum(clean_all), num_updates);
    fprintf('===========================================================================\n');
    fprintf('Estimate = mean over the switch-free windows. Max difference: over those windows / over all windows (switch transients included).\n');
    fprintf('Quantity |        Actual         |       Estimate        | Accuracy | Max difference [Ohm]\n');
    fprintf('R_n      | %12.4f Ohm       | %12.4f Ohm       | %6.1f %% | %.4f / %.4f\n', ...
            Rn_act, Rn_est, acc_Rn, max(abs(real(Zn_err(clean_all)))), max(abs(real(Zn_err))));
    fprintf('L_n      | %14.4e H        | %14.4e H        | %6.1f %% | %.4f / %.4f (reactance)\n', ...
            Ln_act, Ln_est, acc_Ln, max(abs(imag(Zn_err(clean_all)))), max(abs(imag(Zn_err))));
    fprintf('Zn       | %8.4f %+8.4fj Ohm | %8.4f %+8.4fj Ohm | %6.1f %% | %.4f / %.4f\n', ...
            real(Zn_nom), imag(Zn_nom), real(Zn_fin), imag(Zn_fin), acc_Zn, max(abs(Zn_err(clean_all))), max(abs(Zn_err)));
    fprintf('(Zn at %g Hz; the max differences compare the reactance at the frequency of each window.)\n', f_nom);
else
    fprintf('  NEUTRAL WIRE IMPEDANCE Zn  (no neutral sensors: In = I1+I2+I3, Zn estimated from the three phases)\n');
    fprintf('===========================================================================\n');
    if strcmp(nfit.mode, 'none')
        fprintf('WARNING: Zn could not be estimated (%s).\n', nfit.reason);
        fprintf('         Zn was taken as 0, so the Zth values above still contain the neutral path (log Vn and In, or switch loads in at least two phases).\n');
    else
        fprintf('Source-consistency fit (%s): %d switch-free windows, %d iterations, residual %.2f %% of the source EMF\n', ...
                nfit.mode, nfit.n, nfit.iters, 100 * nfit.rel_resid);
        fprintf('Spread of the fit (std of R_n / X_n): %.2g / %.2g Ohm\n', nfit.sd_Zn(1), nfit.sd_Zn(2));
        if ~(nfit.rel_resid <= 0.015)
            fprintf('WARNING: residual above 1.5 %%: the sources do not behave like constant-amplitude sources here (or a switch was missed). Treat Zn with care.\n');
        end
    end
    fprintf('Accuracy = 100 * (1 - |estimate - actual| / |actual|). Max difference = largest |estimate - actual| over the frequency ramp.\n');
    fprintf('Quantity |        Actual         |       Estimate        | Accuracy | Max difference [Ohm]\n');
    fprintf('R_n      | %12.4f Ohm       | %12.4f Ohm       | %6.1f %% | %.4f\n', Rn_act, Rn_est, acc_Rn, max(abs(real(Zn_err))));
    fprintf('L_n      | %14.4e H        | %14.4e H        | %6.1f %% | %.4f (reactance)\n', Ln_act, Ln_est, acc_Ln, max(abs(imag(Zn_err))));
    fprintf('Zn       | %8.4f %+8.4fj Ohm | %8.4f %+8.4fj Ohm | %6.1f %% | %.4f\n', ...
            real(Zn_nom), imag(Zn_nom), real(Zn_fin), imag(Zn_fin), acc_Zn, max(abs(Zn_err)));
    fprintf('(Zn at %g Hz.)\n', f_nom);
end
fprintf('===========================================================================\n');

%% FIGURES
% Continuous Zth graphs: start at the start value; grey bands mark the switch transients (line interpolated there)
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

% Thevenin inductance in henries: L_th = X_th/(2*pi*f) at the frequency of each measurement
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

% Detection view: per-update change of R and X against the thresholds
figure('Name', 'Switch detection from Z_load, 3 phases', 'NumberTitle', 'off', 'Position', [120 40 1150 820]);
for ph = 1:n_ph
    r = res{ph};

    subplot(n_ph, 2, 2*ph - 1);
    plot(t_update(1:end-1), dR_step{ph}, 'b-', 'DisplayName', '\Delta R per update');
    hold on;
    plot([0 T], [ dR_thresh(ph)  dR_thresh(ph)], 'k--', 'DisplayName', 'threshold');
    plot([0 T], [-dR_thresh(ph) -dR_thresh(ph)], 'k--', 'HandleVisibility', 'off');
    plot(r.t_event, zeros(size(r.t_event)), 'rv', 'MarkerFaceColor', 'r', 'DisplayName', 'detected switch');
    plot(r.t_foreign, zeros(size(r.t_foreign)), '^', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], ...
         'DisplayName', 'switch in another phase');
    hold off;
    grid on;
    ylabel('\Delta R (\Omega)');
    title(sprintf('Phase %d: change of R between consecutive windows', ph));
    if ph == 1, legend('Location', 'best', 'FontSize', 7); end
    if ph == n_ph, xlabel('Time (s)'); end
    xlim([0, T]);

    subplot(n_ph, 2, 2*ph);
    plot(t_update(1:end-1), dX_step{ph}, 'r-', 'DisplayName', '\Delta X per update');
    hold on;
    plot([0 T], [ dX_thresh(ph)  dX_thresh(ph)], 'k--', 'DisplayName', 'threshold');
    plot([0 T], [-dX_thresh(ph) -dX_thresh(ph)], 'k--', 'HandleVisibility', 'off');
    plot(r.t_event, zeros(size(r.t_event)), 'rv', 'MarkerFaceColor', 'r', 'DisplayName', 'detected switch');
    plot(r.t_foreign, zeros(size(r.t_foreign)), '^', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], ...
         'DisplayName', 'switch in another phase');
    hold off;
    grid on;
    ylabel('\Delta X (\Omega)');
    title(sprintf('Phase %d: change of X between consecutive windows', ph));
    if ph == 1, legend('Location', 'best', 'FontSize', 7); end
    if ph == n_ph, xlabel('Time (s)'); end
    xlim([0, T]);
end

% Neutral wire impedance: R_n and L_n (measured in every window, or the estimate from the three phases)
figure('Name', 'Neutral wire impedance Zn', 'NumberTitle', 'off', 'Position', [130 70 900 620]);
nz_name  = {'R_n', 'L_n'};
nz_ylab  = {'Resistance (\Omega)', 'Inductance (H)'};
nz_trace = {Rn_calc, Ln_calc};                      % one value per window
nz_est   = [Rn_est, Ln_est];                        % estimate
nz_act   = [Rn_act, Ln_act];                        % actual values (comparison only)
nz_col   = {'b', 'm'};
for p = 1:2
    subplot(2, 1, p);
    hold on;
    if neutral_sensors
        plot(t_update, nz_trace{p}, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0, 'DisplayName', 'every window');
        plot(t_update(clean_all), nz_trace{p}(clean_all), [nz_col{p} '.'], 'MarkerSize', 8, 'DisplayName', 'switch-free windows');
        vals = [nz_trace{p}(clean_all); nz_est(p); nz_act(p)];  % y-range from the switch-free values (transient spikes may leave the plot)
        plot([0 T], nz_est(p) * [1 1], [nz_col{p} '-'], 'LineWidth', 2.5, ...
             'DisplayName', ['Estimated ' nz_name{p} sprintf(' (%.4g), mean of the switch-free windows', nz_est(p))]);
        title(['Neutral wire ' nz_name{p} ' = V_n / I_n']);
    else
        vals = [nz_est(p); nz_act(p)];
        plot([0 T], nz_est(p) * [1 1], [nz_col{p} '-'], 'LineWidth', 2.5, ...
             'DisplayName', ['Estimated ' nz_name{p} sprintf(' (%.4g), from the three phases', nz_est(p))]);
        title(['Neutral wire ' nz_name{p} ' estimated without neutral sensors']);
    end
    plot([0 T], nz_act(p) * [1 1], 'k--', 'LineWidth', 2.0, 'DisplayName', ['Actual ' nz_name{p} sprintf(' (%.4g)', nz_act(p))]);
    span = max(max(vals) - min(vals), 0.5 * nz_act(p));
    ylim(0.5 * (max(vals) + min(vals)) + 0.65 * span * [-1 1]);
    hold off;
    grid on;
    ylabel(nz_ylab{p});
    legend('Location', 'best', 'FontSize', 8);
    if p == 2, xlabel('Time (s)'); end
    xlim([0, T]);
end

%% LOCAL FUNCTIONS

function Zb = branch_z(R, L, C, w)
    % Impedance of a series R-L-C branch at angular frequency w (L = 0: no inductor, C = Inf: no capacitor).
    Zb = R + 1i * w * L;
    if isfinite(C), Zb = Zb - 1i / (w * C); end
end

function Zb = branch_z_solver(R, L, C, w, h, scheme)
    % The same branch as a fixed-step local solver (step h) implements it. The solver replaces the derivative (s = j*w) by a
    % difference formula; with z1 = exp(-j*w*h):
    %   'BE' (Backward Euler)    s = (1 - z1)/h                ~ j*w + w^2*h/2     -> inductor: + w^2*h*L/2 in series;
    %                                                                                   capacitor: + h/(2*C) in series
    %   'TR' (Trapezoidal rule)  s = (2/h)*(1 - z1)/(1 + z1)   = (2/h)*j*tan(w*h/2) -> purely reactive, no damping
    %   anything else            s = j*w (exact)
    z1 = exp(-1i * w * h);
    if strcmpi(scheme, 'BE')
        s = (1 - z1) / h;
    elseif strcmpi(scheme, 'TR')
        s = (2 / h) * (1 - z1) / (1 + z1);
    else
        s = 1i * w;
    end
    Zb = R + s * L;
    if isfinite(C), Zb = Zb + 1 / (s * C); end
end

function [ev_first, ev_last, dR_step, dX_step] = detect_switches(Z_calc, cfg)
    % Switch detection from the changes in Z_load.
    % Z_calc: Z_load of one phase (num_updates x 1). Returns the first and last window touched by each switch (column vectors)
    % and the change of R and X between consecutive windows (for the detection figure). The real and imaginary parts are
    % tested separately (not the magnitude of the complex change); an increase or a decrease of either part counts.
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
end

function r = thevenin_from_load_changes(Vk_load, Vnk_vec, Ik_vec, Nw_vec, tc_vec, f_est, Tw_vec, t_update, ...
                                        ev_first, ev_last, fo_first, fo_last, Zth_true, Zth_start, cfg)
    % Thevenin impedance of one phase, with the neutral drop removed.
    % Inputs: the coefficient sums of one phase (Vk_load: voltage across the load, Ik_vec: current) and of the voltage across
    % the neutral impedance (Vnk_vec), all num_updates x 1, and the window quantities shared by the phases.
    % ev_first / ev_last: windows touched by the switches of this phase. fo_first / fo_last: windows touched by the switches
    % of the other phases (they reach this phase through the neutral). Zth_true (NaN when unknown) and Zth_start belong to
    % this phase. cfg holds the settings.
    num_updates = numel(t_update);
    T = cfg.T;

    % 1. Switching instants of this phase and of the other two
    n_ev      = numel(ev_first);
    t_event   = 0.5 * (t_update(ev_first) + Tw_vec(ev_first) + t_update(ev_last));   % estimated switching instants (s)
    t_foreign = 0.5 * (t_update(fo_first) + Tw_vec(fo_first) + t_update(fo_last));   % the same for the other phases

    % 2. Windows untouched by a switch, and the load state of each window
    is_clean = true(num_updates, 1);
    t_gap    = zeros(n_ev, 2);                        % [start end] (s) of the unclean gap around each switch of this phase
    for e = 1:n_ev
        t_gap(e, :) = [t_update(max(ev_first(e) - cfg.n_guard, 1)), t_update(ev_last(e)) + cfg.t_settle];
        is_clean(t_update >= t_gap(e, 1) & t_update < t_gap(e, 2)) = false;
    end
    for j = 1:numel(fo_first)                         % switches of the other phases: shorter gap, they act only through Vn
        g0 = t_update(max(fo_first(j) - cfg.n_guard, 1));
        g1 = t_update(fo_last(j)) + cfg.t_settle_f;
        is_clean(t_update >= g0 & t_update < g1) = false;
    end
    state = 1 + sum((1:num_updates)' > ev_last', 2);  % 1 = before the first switch of this phase, 2 = after it, ...

    % 3. Phasors from the coefficient sums, source rotation removed
    % The sums refer to t = 0 and grow with the window length. Both cancel in the ratio V/I of one window but not in
    % P = I1/I2 (two different windows), so every window is referred to its own centre and normalised by its length first.
    rot = exp(1i * 2 * pi * f_est .* tc_vec);
    Iph = (Ik_vec ./ Nw_vec) .* rot;                  % current of this phase
    Vlp = (Vk_load ./ Nw_vec) .* rot;                 % voltage across the load
    Vph = ((Vk_load + Vnk_vec) ./ Nw_vec) .* rot;     % load voltage plus the drop in the neutral = voltage at the source terminals
    Eph = Vph + Zth_start * Iph;                      % proxy of the source EMF E = V + Zth*I (with the start value of Zth)

    % The two states of a switch must be compared at the same source phase, whatever the grid frequency did in between.
    % The source phase is read from the phase of the EMF proxy Eph in the clean stretches and carried across the unclean
    % gap at each switch by a smoothing spline; no frequency model is assumed. The proxy is used instead of the phase of the
    % current because a load step moves the phase of the current but not that of the EMF. Zth_start only has to be roughly
    % right (its error enters as dZth*dI/E, below 1 mrad).
    %
    % 3a. Phase of Eph, unwrapped inside every clean stretch. f_est supplies only the whole-cycle count between consecutive
    %     windows (it sits on a 0.1 Hz grid, so integrating it directly would drift by several hundredths of a radian).
    %     Each stretch is anchored to a linear reference phase (ph_ref) so that the numbers stay small.
    ph_ref = 2 * pi * mean(f_est) * (tc_vec - tc_vec(1));
    U = nan(num_updates, 1);
    for s = 1:n_ev+1
        kk = find(is_clean & state == s);
        if numel(kk) < 2, continue; end
        dtc  = diff(tc_vec(kk));
        dphi = angle(Eph(kk(2:end)) ./ Eph(kk(1:end-1)));
        ncyc = round(0.5 * (f_est(kk(2:end)) + f_est(kk(1:end-1))) .* dtc - dphi / (2*pi));
        U(kk) = ph_ref(kk(1)) + angle(Eph(kk(1)) * exp(-1i * ph_ref(kk(1)))) + [0; cumsum(dphi + 2*pi*ncyc)];
    end
    ok   = ~isnan(U);
    rows = find(ok);
    nr   = numel(rows);

    % 3b. Smoothing spline with one free offset per load state:
    %         minimise  sum_clean ( theta_k + c_state(k) - U_k )^2  +  lambda * sum ( 3rd divided difference of theta )^2
    %     The offsets take up any phase step between two load states, so theta follows the source only and is one continuous
    %     curve through every switching boundary. The third-derivative penalty leaves constant frequency and linear ramps
    %     unpenalised; random or non-linear wander slower than f_smooth is followed. Across each unclean gap theta is
    %     interpolated with minimum third-derivative energy. lambda = 36*(2*pi*f_smooth)^-6 puts the half-power point of
    %     the smoother at f_smooth (36 = (3!)^2).
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
        s_off(s_off == s_ref) = [];                          % states with a free offset (none if only one state has clean windows)
        M   = [sparse((1:nr)', rows, 1, nr, n), sparse(double(state(rows) == s_off(:)'))];   % [theta; offsets] -> U
        H   = M' * M;
        H(1:n, 1:n) = H(1:n, 1:n) + lambda * (Dd' * Dd);
        rhs = M' * (U(rows) - ph_ref(rows));
        z   = H \ rhs;
        z   = z + H \ (rhs - H * z);                         % one refinement step (the system is badly scaled)
        theta_src    = z(1:n) + ph_ref;                      % continuous source phase (up to a constant)
        fit_rms_mrad = 1e3 * sqrt(mean((M * z - (U(rows) - ph_ref(rows))).^2));
    end
    Idr = Iph .* exp(-1i * theta_src);                       % current with the source rotation removed
    Vdr = Vph .* exp(-1i * theta_src);                       % voltage at the source terminals, likewise
    Vld = Vlp .* exp(-1i * theta_src);                       % voltage across the load, likewise

    % Grid frequency of every window: the slope of the source phase. f_est cannot be used: it sits on a 0.1 Hz grid, while
    % the two states of a switch differ by only a few tenths of a hertz (see zth_two_states).
    f_src = zeros(num_updates, 1);
    f_src(2:end-1) = (theta_src(3:end) - theta_src(1:end-2)) ./ (tc_vec(3:end) - tc_vec(1:end-2)) / (2 * pi);
    f_src(1)   = f_src(2);
    f_src(end) = f_src(end-1);

    % 4. Apply the formula at every switch of this phase
    Zth_event      = complex(nan(n_ev, 1), nan(n_ev, 1));
    Zth_blind      = complex(nan(n_ev, 1), nan(n_ev, 1));   % the same switch without the neutral correction (for the report)
    Zth_plain      = complex(nan(n_ev, 1), nan(n_ev, 1));   % the plain formula, Zth taken equal in both states (for the report)
    df_event       = nan(n_ev, 1);   % grid frequency after minus before the switch (Hz)
    P_event        = complex(nan(n_ev, 1), nan(n_ev, 1));
    Zth_true_event = complex(nan(n_ev, 1), nan(n_ev, 1));
    sig_R_event    = nan(n_ev, 1);   % predicted error of Re(Zth) from each switch
    sig_X_event    = nan(n_ev, 1);   % predicted error of Im(Zth) from each switch
    f_event        = nan(n_ev, 1);   % frequency at which each switch measured Zth (Hz)
    n_pre          = zeros(n_ev, 1); % clean windows found just before each switch
    n_post         = zeros(n_ev, 1); % clean windows found after its transient

    for e = 1:n_ev
        k_pre  = find(is_clean & state == e,   cfg.n_avg, 'last');    % previous state, just before the switch
        k_post = find(is_clean & state == e+1, cfg.n_avg, 'first');   % current state, after the transient
        n_pre(e)  = numel(k_pre);
        n_post(e) = numel(k_post);
        if n_pre(e) < cfg.n_min || n_post(e) < cfg.n_min, continue; end   % too few clean windows (e.g. start-up)

        I1  = mean(Idr(k_post));                      % current after the switch  (state 1 of the formula)
        I2  = mean(Idr(k_pre));                       % current before the switch (state 2)
        ZL1 = mean(Vld(k_post)) / I1;                 % load alone, after / before the switch
        ZL2 = mean(Vld(k_pre))  / I2;
        if abs(real(ZL1 - ZL2)) <= cfg.dR_thresh && abs(imag(ZL1 - ZL2)) <= cfg.dX_thresh
            continue;                                 % no net load change: false alarm
        end

        Zb1 = mean(Vdr(k_post)) / I1;                 % load + neutral drop, as seen from the source terminals
        Zb2 = mean(Vdr(k_pre))  / I2;
        P   = I1 / I2;                                % P = I_1 / I_2
        f_post = mean(f_src(k_post));                 % grid frequency of the two states (the grid moves between them)
        f_pre  = mean(f_src(k_pre));
        P_event(e)        = P;
        Zth_event(e)      = zth_two_states(I1, I2, Zb1, Zb2, f_post, f_pre);   % Zth = (Z_L2 - P*Z_L1)/(P - 1), Z_L = (V_load + Vn)/I
        Zth_blind(e)      = zth_two_states(I1, I2, ZL1, ZL2, f_post, f_pre);   % the same with Z_L = V_load/I: includes the neutral path
        Zth_plain(e)      = (Zb2 - P * Zb1) / (P - 1);                         % the formula as it stands, Zth equal in both states
        df_event(e)       = f_post - f_pre;
        Zth_true_event(e) = mean(Zth_true(k_post));
        f_event(e)        = f_post;

        % error gain: dZth = err_gain * (amplitude error + j*phase error of P)
        err_gain = (Zb1 - Zb2) * P / (P - 1)^2;
        sig_R_event(e) = hypot(cfg.sig_amp * real(err_gain), cfg.sig_phase * imag(err_gain));
        sig_X_event(e) = hypot(cfg.sig_amp * imag(err_gain), cfg.sig_phase * real(err_gain));
    end

    % 4b. Update only the part of Zth that each switch measured reliably
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

    % Inductance view of the reactance: L = X/(2*pi*f) at the frequency of each measurement (L does not change with the
    % frequency ramp, X does). Same switches, weights and start value as X_th.
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

    r = struct('n_ev', n_ev, 't_event', t_event, 't_foreign', t_foreign, 't_gap', t_gap, 'is_clean', is_clean, ...
               'state', state, 'fit_rms_mrad', fit_rms_mrad, 'n_pre', n_pre, 'n_post', n_post, ...
               'Zth_event', Zth_event, 'Zth_blind', Zth_blind, 'Zth_plain', Zth_plain, 'df_event', df_event, ...
               'P_event', P_event, 'Zth_true_event', Zth_true_event, ...
               'sig_R_event', sig_R_event, 'sig_X_event', sig_X_event, 'use_R', use_R, 'use_X', use_X, ...
               'Zth_final', Zth_final, 'sig_R_final', sig_R_final, 'sig_X_final', sig_X_final, ...
               'Zth_true_final', Zth_true_final, 't_zth', t_zth, 'Zth_calc', Zth_calc, ...
               'f_event', f_event, 'L_event', L_event, 'sig_L_event', sig_L_event, 'L_final', L_final, ...
               'sig_L_final', sig_L_final, 'L_trace', L_trace);
end

function Z = zth_two_states(I1, I2, ZL1, ZL2, f1, f2)
    % Thevenin formula when the grid frequency moves between the two states.
    % I1, ZL1, f1: current, load impedance and grid frequency after the switch (state 1); I2, ZL2, f2: before it (state 2).
    % P = I1/I2, Zth = (ZL2 - P*ZL1)/(P - 1) follows from E = I1*(ZL1 + Zth) = I2*(ZL2 + Zth) with the same Zth in both states.
    % On a frequency ramp the reactance of Zth differs between the states: Zth = R + j*2*pi*f*L with R and L constant. The grid
    % moves a few tenths of a hertz between the states, which changes the reactance of a 0.2 Ohm source by about 1 mOhm; the
    % small |P - 1| of a weak switch multiplies that by 1/|P - 1| (5 to 10 times) and it lands in R_th. The same equation is
    % therefore solved for R and L:
    %       R*(I1 - I2) + j*L*(w1*I1 - w2*I2) = I2*ZL2 - I1*ZL1       (w = 2*pi*f; a complex equation, two real unknowns)
    % With f1 = f2 this is exactly Zth = (ZL2 - P*ZL1)/(P - 1). The result is Zth at the frequency f1 of the state after the switch.
    w1 = 2 * pi * f1;
    w2 = 2 * pi * f2;
    a  = I1 - I2;
    b  = 1i * (w1 * I1 - w2 * I2);
    c  = I2 * ZL2 - I1 * ZL1;
    RL = [real(a) real(b); imag(a) imag(b)] \ [real(c); imag(c)];
    Z  = RL(1) + 1i * w1 * RL(2);
end

function fit = fit_neutral_from_phases(Vph, Iph, f_k, st, mode_in)
    % Neutral impedance from the three phases (used when Vn and In are not logged).
    % Vph, Iph: n x 3 window phasors (voltage across the load, current) of the switch-free windows of the three phases.
    % f_k: n x 1 frequency of each window. st: n x 3 load state of each phase in each window (1, 2, ... counting its switches).
    % mode_in: 'auto', 'rigid' or 'amplitude' (see rigid_sources at the top).
    %
    % Model. In = I1 + I2 + I3 flows in the neutral, so in every window the source EMF of phase p is
    %       C_p = V_p + Zth_p*I_p + Zn*In ,     Zth_p = R_p + j*X_p*(f/50) ,     Zn = Rn + j*Xn*(f/50) .
    % C_p must not depend on the load. Its size is constant, |C_p| = |E_p|, and for sources that rotate together it is
    % E_p times one unit phasor u that is the same for the three phases (and different in every window: the frequency
    % moves). Unknowns: R_p, X_p, Rn, Xn and the EMFs E_p. For given parameters the best u follows in closed form
    % (u = S/|S|, S = sum_p conj(E_p)*C_p), which leaves a small nonlinear least-squares problem for the parameters,
    % solved by Levenberg-Marquardt. A switch in one phase changes In and so moves C_p of the other phases by Zn*dIn unless
    % Zn is right: that is what fixes Zn. Nothing depends on the grid frequency or on how it varies.
    %
    % Zn can only be told apart from the source impedances if the loads change in at least two phases inside the windows
    % given here and the fit is well determined; otherwise the result is mode 'none' (Zn = 0) with the reason in fit.reason
    % (a wrong Zn would be worse than none).
    n = size(Vph, 1);
    fit = struct('Rn', 0, 'Xn', 0, 'R', zeros(1, 3), 'X', zeros(1, 3), 'E', zeros(1, 3), 'mode', 'none', 'label', 'none', ...
                 'n', n, 'iters', 0, 'rel_resid', NaN, 'rel_resid_rigid', NaN, 'sd_Zn', [NaN NaN], 'reason', '');
    if n < 30
        fit.reason = sprintf('only %d switch-free windows, at least 30 are needed', n);
    elseif sum(max(st, [], 1) > min(st, [], 1)) < 2
        fit.reason = 'the loads of fewer than two phases change inside the switch-free windows, so Zn cannot be told apart from the source impedances';
    end
    if ~isempty(fit.reason)
        warning('Zn cannot be estimated: %s. Zn is taken as 0.', fit.reason);
        return;
    end
    In = sum(Iph, 2);                                 % neutral current
    fr = f_k / 50;                                    % frequency relative to 50 Hz (the reactances scale with it)
    C0 = Vph(1, :);                                   % EMFs of the first window seen with Zth = 0 and Zn = 0 (starting values)
    Escale = mean(abs(C0));                           % the EMF parameters are kept near 1
    if strcmpi(mode_in, 'rigid')
        modes = {'rigid'};
    elseif strcmpi(mode_in, 'amplitude')
        modes = {'amplitude'};
    else
        modes = {'rigid', 'amplitude'};               % 'auto': the second one only if the first fits badly
    end
    for im = 1:numel(modes)
        rigid = strcmp(modes{im}, 'rigid');
        if rigid
            E0 = C0 * exp(-1i * angle(C0(1)));        % common rotation removed: E_1 real
            th = zeros(1, 13);                        % [R1 R2 R3  X1 X2 X3  Rn Xn  EMF parameters]
            th(9:13) = [real(E0(1)), real(E0(2)), imag(E0(2)), real(E0(3)), imag(E0(3))] / Escale;
        else
            th = zeros(1, 11);
            th(9:11) = abs(C0) / Escale;              % amplitudes of the three sources
        end
        resid = @(p) neutral_resid(p, Vph, Iph, In, fr, Escale, rigid);
        r      = resid(th);
        cost   = r' * r;
        lambda = 1e-2;
        n_it   = 0;
        for it = 1:200
            J     = numeric_jacobian(resid, th, r);
            A     = J' * J;
            gvec  = J' * r;
            dA    = diag(A) + 1e-6 * max(diag(A));      % damping scaled to the problem (a parameter that hardly acts must not run away)
            improved = false;
            step  = zeros(numel(th), 1);
            for tries = 1:20
                step = -(A + lambda * diag(dA)) \ gvec;
                tn   = th + step.';
                rn   = resid(tn);
                cn   = rn' * rn;
                if cn < cost
                    th = tn;  r = rn;  cost = cn;
                    lambda = max(lambda / 3, 1e-9);
                    improved = true;
                    break;
                end
                lambda = lambda * 4;
            end
            n_it = it;
            if ~improved || norm(step) < 1e-10, break; end
        end
        if rigid
            E = Escale * [th(9), th(10) + 1i * th(11), th(12) + 1i * th(13)];
        else
            E = Escale * th(9:11);
        end
        % residual relative to the size of the sources: rms of |C_p - E_p*u| (rigid; real and imaginary part are two entries of
        % r, hence the factor 2) or of ||C_p| - |E_p|| (amplitude)
        if rigid
            rel_resid = sqrt(2 * cost / numel(r)) / mean(abs(E));
        else
            rel_resid = sqrt(cost / numel(r)) / mean(abs(E));
        end
        % spread of the fitted parameters from J'*J at the solution (scaled by the column sizes, so that a parameter on which
        % the residual hardly depends shows up as a huge spread instead of being hidden by the conditioning of the matrix)
        J     = numeric_jacobian(resid, th, r);
        A     = J' * J;
        s_A   = sqrt(diag(A));
        sd_th = sqrt(diag(inv(A ./ (s_A * s_A')))) ./ s_A * sqrt(cost / max(numel(r) - numel(th), 1));
        fit.Rn = th(7);                               % [Ohm]
        fit.Xn = th(8);                               % [Ohm] at 50 Hz
        fit.R  = th(1:3);                             % source resistances found on the way (by-product, not used)
        fit.X  = th(4:6);                             % source reactances at 50 Hz (by-product, not used)
        fit.E  = E;                                   % source EMFs (first window's rotation)
        fit.n  = n;
        fit.iters = n_it;
        fit.rel_resid = rel_resid;
        fit.sd_Zn = [sd_th(7), sd_th(8)];             % [Ohm]
        if rigid
            fit.mode = 'rigid';       fit.label = 'sources rotate together';
            fit.rel_resid_rigid = rel_resid;
            if rel_resid <= 0.007 || numel(modes) == 1, break; end   % 'auto': keep it unless the residual exceeds 0.7 %
        else
            fit.mode = 'amplitude';   fit.label = 'amplitude of each source only';
        end
    end

    % Well determined? The spread of Rn and Xn must be small against the load impedances (the scale of the problem).
    % Written so that a NaN spread (singular J'*J) also fails.
    scale_Z = median(abs(Vph(:) ./ Iph(:)));
    if ~(max(fit.sd_Zn) <= 0.02 * scale_Z)
        fit.reason = sprintf('the fit is poorly determined (std of R_n / X_n: %.3g / %.3g Ohm against loads of about %.3g Ohm)', ...
                             fit.sd_Zn(1), fit.sd_Zn(2), scale_Z);
        warning('Zn cannot be estimated: %s. Zn is taken as 0.', fit.reason);
        fit.Rn = 0;
        fit.Xn = 0;
        fit.mode = 'none';
        fit.label = 'none';
    end
end

function J = numeric_jacobian(resid, th, r)
    % Numerical Jacobian of a residual function (forward differences).
    % resid: function handle returning the residual column vector for a parameter row vector; r = resid(th).
    J = zeros(numel(r), numel(th));
    for j = 1:numel(th)
        h  = 1e-6 * max(1, abs(th(j)));
        tp = th;  tp(j) = tp(j) + h;
        J(:, j) = (resid(tp) - r) / h;
    end
end

function r = neutral_resid(th, V, I, In, fr, Escale, rigid)
    % Residual of the neutral fit (see fit_neutral_from_phases).
    % th = [R1 R2 R3  X1 X2 X3  Rn Xn  EMF parameters]; V, I: n x 3; In, fr: n x 1. Returns a real column vector.
    Zp = th(1:3) + 1i * th(4:6) .* fr;                % n x 3 source impedances at the frequency of each window
    Zn = th(7) + 1i * th(8) * fr;                     % n x 1 neutral impedance
    C  = V + Zp .* I + Zn .* In;                      % n x 3: source EMF of each phase, as seen in every window
    if rigid
        E   = Escale * [th(9), th(10) + 1i * th(11), th(12) + 1i * th(13)];   % 1 x 3
        S   = sum(conj(E) .* C, 2);                   % the best common rotation is u = S/|S|
        res = C - E .* (S ./ abs(S));
        r   = [real(res(:)); imag(res(:))];
    else
        E   = Escale * th(9:11);                      % 1 x 3 amplitudes
        res = abs(C) - abs(E);
        r   = res(:);
    end
end

function [A6, D6, D5] = dwt_sym4_l6(sig, LoD, HiD)
    % 6-level sym4 DWT with periodic extension. Returns the level-6 approximation A6 and the level-6 and level-5 details D6, D5.
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

        for k = 0:7                                  % filter sums for all n at once
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