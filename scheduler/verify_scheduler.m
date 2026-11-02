%% VERIFY_SCHEDULER.m
% Verification: phase scheduler fires correctly with hover disabled.
%
% To test WITH hover: set HOVER_ENABLED = true below AND set
%                     params.n_phases = 4 (done automatically below).
%

clear; clc;

script_dir   = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(genpath(project_root));

HOVER_ENABLED = false;   % set true (and params.n_phases=4 below) to run
%                          the 4-phase hover test instead

%% Params struct 
params          = struct();
params.n_phases = 3 + HOVER_ENABLED;   % 3 when hover off, 4 when hover on
% Vehicle mass for synthetic trajectory (not used by scheduler logic,
% but required so x_synth(17,:) is physically reasonable)
params.m0 = 16.88;   % kg — total initial mass (dry + propellant)
params.ms = 8.00;    % kg — dry mass

%% Nominal phase times
if HOVER_ENABLED
    t_phase_nominal = [12; 17; 30; 35];   % 4-element: [t1 t2 t3 t4] end times
else
    t_phase_nominal = [12; 28; 33];       % 3-element: [t1 t3 t4] end times
end

%% Synthetic 17-state trajectory
% State layout: [r(3); v(3); q(4); b_g(3); b_a(3); m(1)]  (17x1)
%   Scheduler only reads r(1:3) and v(4:6); all other states can be nominal.
dt_sched = 0.5;
T_sim    = 36;
t_vec    = 0 : dt_sched : T_sim;
N        = numel(t_vec);

pos_start = [0;   0;    0  ];
pos_hover = [100; -100; 300];
pos_desc  = [200; -200; 5  ];
pos_land  = [200; -200; 0  ];

if HOVER_ENABLED
    t1e = t_phase_nominal(1);
    t2e = t_phase_nominal(2);
    t3e = t_phase_nominal(3);
    t4e = t_phase_nominal(4);
else
    t1e = t_phase_nominal(1);
    t2e = t1e;       % no hover (zero hover duration)
    t3e = t_phase_nominal(2);
    t4e = t_phase_nominal(3);
end

x_synth = zeros(17, N);   % 17-state EKF vector
for k = 1:N
    tk = t_vec(k);

    % --- Position and velocity (scheduler-relevant states) ---
    if tk <= t1e
        alpha = tk / t1e;
        r = pos_start + alpha * (pos_hover - pos_start);
        v = (pos_hover - pos_start) / t1e;
        if alpha > 0.8
            v = v * (1 - (alpha-0.8)/0.2 * 0.95);   % decelerate near apex
        end
    elseif HOVER_ENABLED && tk <= t2e
        r = pos_hover;
        v = zeros(3,1);
    elseif tk <= t3e
        alpha = (tk - t2e) / (t3e - t2e);
        r = pos_hover + alpha * (pos_desc - pos_hover);
        v = (pos_desc - pos_hover) / (t3e - t2e);
        if alpha > 0.85
            v = v * (1 - (alpha-0.85)/0.15 * 0.9);  % decelerate near landing
        end
    elseif tk <= t4e
        alpha = (tk - t3e) / (t4e - t3e);
        r = pos_desc + alpha * (pos_land - pos_desc);
        v = (pos_land - pos_desc) / (t4e - t3e);
    else
        r = pos_land;
        v = zeros(3,1);
    end

    x_synth(1:3, k) = r;
    x_synth(4:6, k) = v;

    % --- Attitude: upright throughout (quaternion [1;0;0;0]) ---
    x_synth(7,  k) = 1;   % qw=1

    % --- Biases: zero (nominal) ---
    % x_synth(11:16, k) = 0;  (already zero from zeros initialisation)

    % --- Mass: linear depletion m0 → ~0.52*m0 over mission ---
    % Approximate; exact value doesn't affect scheduler logic.
    t_total = t4e;
    x_synth(17, k) = params.m0 * (1 - 0.48 * min(tk, t_total) / t_total);
end

%% Run nominal scheduler
sched = sched_state_init(0);
phase_log = zeros(1, N);
t_in_log  = zeros(1, N);
t_to_log  = zeros(1, N);
trans_times = struct('ph1', NaN, 'ph2', NaN, 'ph3', NaN, 'ph4', NaN);

for k = 1:N
    [pid, t_in, t_to, sched] = phase_scheduler(t_vec(k), x_synth(:,k), ...
                                                t_phase_nominal, params, sched);
    phase_log(k) = pid;
    t_in_log(k)  = t_in;
    t_to_log(k)  = t_to;

    switch pid
        case 1, if isnan(trans_times.ph1), trans_times.ph1 = t_vec(k); end
        case 2, if isnan(trans_times.ph2), trans_times.ph2 = t_vec(k); end
        case 3, if isnan(trans_times.ph3), trans_times.ph3 = t_vec(k); end
        case 4, if isnan(trans_times.ph4), trans_times.ph4 = t_vec(k); end
    end
end

fprintf('Phase entry times:  1=%.1f  2=%s  3=%.1f  4=%.1f\n', ...
    trans_times.ph1, ...
    ternary(isnan(trans_times.ph2), 'N/A (hover off)', sprintf('%.1f',trans_times.ph2)), ...
    trans_times.ph3, trans_times.ph4);

%% Tests
all_pass = true;

% Test 1: Phase 1-3 timing
fprintf('\n--- Test 1: Phase 1-3 timing ---\n');
err13 = abs(trans_times.ph3 - t1e);
if err13 <= 1.0
    fprintf('  PASS: phase 3 entered at %.2f s (nominal t1_end=%.1f, err=%.2f s)\n', trans_times.ph3, t1e, err13);
else
    fprintf('  FAIL: phase 3 entered at %.2f s (expected ~%.1f, err=%.2f s)\n', trans_times.ph3, t1e, err13);
    all_pass = false;
end

% Test 2: Phase 3-4 timing
fprintf('\n--- Test 2: Phase 3→4 timing ---\n');
err34 = abs(trans_times.ph4 - t3e);
if err34 <= 1.0
    fprintf('  PASS: phase 4 entered at %.2f s (nominal t3_end=%.1f, err=%.2f s)\n', trans_times.ph4, t3e, err34);
else
    fprintf('  FAIL: phase 4 entered at %.2f s (expected ~%.1f, err=%.2f s)\n', trans_times.ph4, t3e, err34);
    all_pass = false;
end

% Test 3: Forward-only
fprintf('\n=== Test 3: Forward-only phase order ===\n');
if all(diff(phase_log) >= 0)
    fprintf('  PASS: phase IDs never decrease\n');
else
    fprintf('  FAIL: phase went backwards\n');
    all_pass = false;
end

% Test 4: Phase 2 never appears when hover disabled
fprintf('\n=== Test 4: Phase 2 absent (hover disabled) ===\n');
if ~HOVER_ENABLED
    if ~any(phase_log == 2)
        fprintf('  PASS: phase 2 (hover) never activated\n');
    else
        fprintf('  FAIL: phase 2 appeared even though hover is disabled\n');
        all_pass = false;
    end
else
    fprintf('  SKIP: hover is enabled — phase 2 expected\n');
end

% Test 5: Hysteresis — state spike before time gate must not trigger
fprintf('\n=== Test 5: Hysteresis spike rejection ===\n');
sched_h = sched_state_init(0);
phase_h = zeros(1, N);
for k = 1:N
    x_k = x_synth(:, k);
    % Inject hover-like position/velocity at t=8 s (before time gate at t1_end-2=10 s)
    if abs(t_vec(k) - 8.0) < 0.25
        x_k(1:3) = pos_hover;
        x_k(4:6) = zeros(3,1);
    end
    [pid, ~, ~, sched_h] = phase_scheduler(t_vec(k), x_k, t_phase_nominal, params, sched_h);
    phase_h(k) = pid;
end
if ~any(phase_h(t_vec < 9) > 1)
    fprintf('  PASS: state spike at t=8 s did not trigger early transition\n');
else
    fprintf('  FAIL: spurious transition before t=9 s\n');
    all_pass = false;
end

% Test 6: t_in_phase resets at transitions
fprintf('\n=== Test 6: t_in_phase resets at transitions ===\n');
trans_idx = find(diff(phase_log) ~= 0) + 1;
pass_reset = true;
for ki = trans_idx
    if t_in_log(ki) > dt_sched + 0.01
        fprintf('  FAIL: t_in_phase = %.3f s at transition (expected ~0)\n', t_in_log(ki));
        pass_reset = false;
        all_pass   = false;
    end
end
if pass_reset
    fprintf('  PASS: t_in_phase resets correctly at all transitions\n');
end

%% Summary
fprintf('\n%s\n', repmat('=',1,45));
if all_pass
    fprintf('  ALL PASS — scheduler ready for Step 4\n');
else
    fprintf('  Some tests FAILED — check transition thresholds\n');
end

%% Plot
figure('Name','Phase Scheduler Verification','Position',[80 80 1100 700]);

subplot(2,3,1);
stairs(t_vec, phase_log, 'b', 'LineWidth', 1.5); hold on;
for ph = [3 4]
    field = sprintf('ph%d', ph);
    if ~isnan(trans_times.(field))
        xline(trans_times.(field), '--r');
        text(trans_times.(field)+0.2, ph-0.15, sprintf('→%d',ph), 'Color','r');
    end
end
xline(t1e, ':k', 'nom t1'); xline(t3e, ':k', 'nom t3');
grid on; xlabel('t (s)'); ylabel('Phase ID');
ylim([0.5 4.5]); yticks([1 2 3 4]);
title(sprintf('Phase ID  (n\\_phases=%d)', params.n_phases));

subplot(2,3,2);
plot(t_vec, t_in_log, 'b'); grid on;
xlabel('t (s)'); ylabel('t_{in\_phase} (s)');
title('Time in current phase');

subplot(2,3,3);
plot(t_vec, x_synth(3,:), 'b'); hold on;
yline(300,'--k','hover alt'); yline(5,':k','desc alt');
grid on; xlabel('t (s)'); ylabel('Z (m)');
title('Altitude profile');

subplot(2,3,4);
plot(t_vec, x_synth(1,:),'b', t_vec, x_synth(2,:),'r', t_vec, x_synth(3,:),'y');
grid on; xlabel('t (s)'); legend('X','Y','Z'); title('Trajectory (XYZ)');

subplot(2,3,5);
stairs(t_vec, phase_h, 'b', 'LineWidth', 1.5); grid on;
xlabel('t (s)'); ylabel('Phase ID');
ylim([0.5 4.5]); yticks([1 2 3 4]);
title('Test 5: phase under spike injection');

subplot(2,3,6);
plot(t_vec, x_synth(17,:), 'b'); grid on;
xlabel('t (s)'); ylabel('Mass (kg)');
title('Synthetic mass profile');

sgtitle(sprintf('Step 3 — Phase Scheduler Verification  (n\\_phases=%d)', params.n_phases));


function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end