function test_modeB_descent
% TEST_MODEB_CLOSED_LOOP  Pure-MATLAB closed-loop test of Mode B landing.
%
% Runs phase 3 (the descent manoeuvre) in closed loop without Simulink:
%   - Initial state taken directly from the reference's phase-3-entry node
%   - 3-DoF translational dynamics integrated by ode45 in MATLAB
%   - Mode B MPC replans every params.dt_guidance seconds
%   - Trajectory buffer between replans
%   - No attitude loop, no MEKF, "perfect" thrust-direction tracking


% L. White, University of Sheffield, May 2026

clear; clc; close all;
addpath(genpath('.'));

fprintf('  TEST_MODEB_DESCENT - Step 9 \n');

t_start_wallclock = tic;

%% Setup
params = plant_params();
ref    = params.reference;

% 2 Hz guidance rate (default for landing)
params.dt_guidance = 0.5;

if ~isfield(params, 'socp_optimizer') || isempty(params.socp_optimizer)
    fprintf('Compiling SOCP (one-time) ...\n');
    params.socp_optimizer = build_socp(params);
end

%% Initial state: Phase-3 entry, sampled from reference
t_phase3_start = ref.phase_transition_times(3);
t_phase3_end   = ref.phase_transition_times(4);   % = Phase 4 entry

[t_unique, ia] = unique(ref.t, 'stable');
x_unique = ref.x(ia, :);

% Initial state
x_at_entry = interp1(t_unique, x_unique, t_phase3_start, 'linear').';
x0_si = [x_at_entry(1:3); x_at_entry(4:6); x_at_entry(7)];

% Target: Phase 4 entry state (position + velocity, mass left free)
x_at_transit = interp1(t_unique, x_unique, t_phase3_end, 'linear').';
target_state = x_at_transit([1:3, 4:6]);   % 6x1: [r_target; v_target]

fprintf('Phase-3 entry conditions:\n');
fprintf('  t_phase3_start = %.3f s\n', t_phase3_start);
fprintf('  r  = [%7.2f %7.2f %7.2f] m\n',   x0_si(1:3));
fprintf('  v  = [%7.3f %7.3f %7.3f] m/s\n', x0_si(4:6));
fprintf('  m  = %.3f kg\n\n',                x0_si(7));
fprintf('Phase-4 entry (target):\n');
fprintf('  r  = [%7.2f %7.2f %7.2f] m\n',   target_state(1:3));
fprintf('  v  = [%7.3f %7.3f %7.3f] m/s\n', target_state(4:6));
fprintf('  reference duration = %.3f s\n\n', t_phase3_end - t_phase3_start);

%% Simulation timing
dt_sim = 0.01;

% Descent is free-final-time. Cap sim length at sigma_high + safety.
sigma_max_s        = params.scp.sigma_desc_high;
t_sim_end          = t_phase3_start + sigma_max_s + 1.0;
t_no_replan_after  = t_sim_end - 1.0;
fprintf('Max simulation duration: %.1f s\n\n', t_sim_end - t_phase3_start);

t_grid = t_phase3_start:dt_sim:t_sim_end;
N_sim  = numel(t_grid);
%% Mode B state
prev_sol        = [];
rejection_count = 0;
fallback_active = false;
last_mpc_t      = -inf;
buf_t = []; buf_T = []; buf_n = [];

%% Logging
x_log     = zeros(7, N_sim);
T_cmd_log = zeros(1, N_sim);
mpc_calls = struct('t',{}, 'valid',{}, 'fault',{}, 'exit_code',{}, ...
                   'iter',{}, 'wall_ms',{}, 'lcvx_N',{}, ...
                   'sigma_sol_s',{}, 'pos_err',{}, 'vel_err',{}, ...
                   'reject',{}, 'fb',{});
n_mpc_calls = 0;
hit_imax = false;

%% Main loop
x = x0_si;

fprintf('Per-call diagnostic:\n');
fprintf('  t[s]   sigma  iter  valid fault exit  lcvx[N]   pos_err  vel_err  rej  fb\n');
fprintf('  %s\n', repmat('-', 1, 76));

for k = 1:N_sim
    t = t_grid(k);

    %% Replan decision
    do_replan = (isempty(buf_t) || (t - last_mpc_t >= params.dt_guidance)) ...
                && (t < t_no_replan_after) ...
                && ~fallback_active;

    if do_replan
        x_hat = zeros(17,1);
        x_hat(1:6)  = x(1:6);
        x_hat(7:10) = [1; 0; 0; 0];
        x_hat(17)   = x(7);

        
        sigma_bounds = [params.scp.sigma_desc_low, params.scp.sigma_desc_high];

        t_call = tic;
        [~, ~, sol_out, rejection_count, fb] = mpc_solve( ...
            x_hat, t, target_state, sigma_bounds, ...
            prev_sol, rejection_count, params);
        t_call_ms = toc(t_call) * 1000;

        last_mpc_t  = t;
        n_mpc_calls = n_mpc_calls + 1;

        % Compute terminal errors the SCP is promising — against the
        % descent target (Phase-4 entry), not the landing pad.
        pos_err_m   = norm(sol_out.x_k(1:3,end)*params.U_L - target_state(1:3));
        vel_err_mps = norm(sol_out.x_k(4:6,end)*params.U_V - target_state(4:6));

        % Pull final-iter exit code from history if available
        exit_code = NaN;
        if isfield(sol_out, 'iter_history') && ~isempty(sol_out.iter_history)
            exit_code = sol_out.iter_history(end).exit_code;
        end

        if sol_out.iter >= params.scp.i_max
            hit_imax = true;
        end

        mpc_calls(end+1) = struct( ...
            't',          t - t_phase3_start, ...
            'valid',      sol_out.valid, ...
            'fault',      sol_out.fault_code, ...
            'exit_code',  exit_code, ...
            'iter',       sol_out.iter, ...
            'wall_ms',    t_call_ms, ...
            'lcvx_N',     sol_out.lcvx_slack_N, ...
            'sigma_sol_s', sol_out.sigma * params.U_T, ...
            'pos_err',    pos_err_m, ...
            'vel_err',    vel_err_mps, ...
            'reject',     rejection_count, ...
            'fb',         fb); %#ok<AGROW>

        fprintf('  %5.2f  %5.2f  %4d  %5d %5d %4d  %7.2e %7.3f  %6.3f  %3d  %2d\n', ...
            t - t_phase3_start, sol_out.sigma * params.U_T, sol_out.iter, sol_out.valid, ...
            sol_out.fault_code, exit_code, sol_out.lcvx_slack_N, ...
            pos_err_m, vel_err_mps, rejection_count, fb);

        if sol_out.valid
            t_splice = sol_out.t_grid_si(1);
            keep_idx = buf_t < t_splice;
            buf_t = [buf_t(keep_idx); sol_out.t_grid_si(:)];
            buf_T = [buf_T(keep_idx); sol_out.T_k_si(:)];
            buf_n = [buf_n(keep_idx, :); sol_out.n_hat_k.'];
            [buf_t, si] = sort(buf_t);
            buf_T = buf_T(si);
            buf_n = buf_n(si, :);
            prev_sol = sol_out;
        end

        if fb && ~fallback_active
            fallback_active = true;
            fprintf('  ** fallback latched at t=%.2f (3 consecutive rejects) **\n', ...
                t - t_phase3_start);
            fprintf('  Continuing — Phase-4 entry detector may still fire shortly.\n');
        end
    end

    %% Read command
    if ~isempty(buf_t) && ~fallback_active && t >= buf_t(1)
        t_clamped = min(buf_t(end), t);
        T_cmd = interp1(buf_t, buf_T, t_clamped, 'linear');
        n_hat = interp1(buf_t, buf_n, t_clamped, 'linear').';
        if norm(n_hat) > 1e-6, n_hat = n_hat / norm(n_hat); else, n_hat = [0;0;1]; end
        T_cmd = max(params.Tmin, min(params.Tmax, T_cmd));
    else
        x_hat_A = zeros(17, 1);
        x_hat_A(1:6)  = x(1:6);
        x_hat_A(7:10) = [1; 0; 0; 0];
        x_hat_A(17)   = x(7);
        [T_cmd, n_hat, ~] = guidance_track_reference(t, x_hat_A, ref, params);
        n_hat = n_hat(:);
        if norm(n_hat) > 1e-6, n_hat = n_hat / norm(n_hat); else, n_hat = [0;0;1]; end
        T_cmd = max(params.Tmin, min(params.Tmax, T_cmd));
    end

    T_vec = T_cmd * n_hat;
    Gamma = T_cmd;

    x_log(:, k)   = x;
    T_cmd_log(k)  = T_cmd;
    
    %% Phase 4 entry detection (build plan §12.10 Test 9a)
    r_err = norm(x(1:3) - target_state(1:3));
    v_err = norm(x(4:6) - target_state(4:6));
    if r_err < 3.0 && v_err < 0.5 && t > t_phase3_start + 2.0
        fprintf('  PHASE 4 ENTRY at t=%.2f s  (r_err=%.2f m, v_err=%.3f m/s)\n', ...
            t - t_phase3_start, r_err, v_err);
        x_log     = x_log(:, 1:k);
        T_cmd_log = T_cmd_log(1:k);
        t_grid    = t_grid(1:k);
        break;
    end

    %% Integrate plant forward
    odefun = @(~, xv) transl_dynamics(xv, T_vec, Gamma, params);
    [~, X] = ode45(odefun, [t, t+dt_sim], x, ...
                   odeset('RelTol', 1e-6, 'AbsTol', 1e-8));
    x = X(end, :).';
end

t_wallclock = toc(t_start_wallclock);
fprintf('\nClosed loop done in %.1f s wall-clock\n', t_wallclock);

%%  DESCENT SUMMARY 
final_pos    = x_log(1:3, end);
final_vel    = x_log(4:6, end);
final_mass   = x_log(7, end);
r_err_final  = norm(final_pos - target_state(1:3));
v_err_final  = norm(final_vel - target_state(4:6));
fuel_burned  = x0_si(7) - final_mass;

% Reference fuel over Phase 3
ref_mass_at_entry  = interp1(t_unique, x_unique(:,7), t_phase3_start, 'linear');
ref_mass_at_transit= interp1(t_unique, x_unique(:,7), t_phase3_end,   'linear');
ref_fuel_burn      = ref_mass_at_entry - ref_mass_at_transit;

valid_calls = [mpc_calls.valid];
if any(valid_calls)
    max_lcvx_N    = max([mpc_calls(valid_calls).lcvx_N]);
    mean_iter     = mean([mpc_calls(valid_calls).iter]);
    solve_times   = [mpc_calls(valid_calls).wall_ms];
    if numel(solve_times) > 1
        max_solve_ms = max(solve_times(2:end));   % exclude first-call YALMIP warmup
    else
        max_solve_ms = solve_times(1);
    end
    sigma_history = [mpc_calls(valid_calls).sigma_sol_s];
    if numel(sigma_history) > 1
        dsigma      = diff(sigma_history);
        sigma_jump  = max(abs(dsigma));               % biggest step
        n_increase  = sum(dsigma > 0.5);              % count of upward jumps > 0.5 s
    else
        sigma_jump  = 0; n_increase = 0;
    end
else
    max_lcvx_N = NaN; mean_iter = NaN; max_solve_ms = NaN;
    sigma_jump = NaN; n_increase = NaN;
    sigma_history = [];
end

fprintf('\n--- DESCENT SUMMARY ---\n');
fprintf('Final r           : [%7.2f %7.2f %7.2f] m\n', final_pos);
fprintf('Target (Ph4 entry): [%7.2f %7.2f %7.2f] m\n', target_state(1:3));
fprintf('Position error    : %.3f m\n', r_err_final);
fprintf('Velocity error    : %.3f m/s\n', v_err_final);
fprintf('Fuel burned       : %.3f kg  (ref: %.3f kg, %+.1f%%)\n', ...
    fuel_burned, ref_fuel_burn, 100*(fuel_burned-ref_fuel_burn)/ref_fuel_burn);
fprintf('Mean SCP iters    : %.2f\n', mean_iter);
fprintf('Max LCvx slack    : %.3e N\n', max_lcvx_N);
fprintf('Max solve time    : %.1f ms\n', max_solve_ms);
fprintf('sigma history (s) : [%s]\n', num2str(sigma_history, '%.2f '));
fprintf('sigma max jump    : %.3f s\n', sigma_jump);
fprintf('sigma up-jumps    : %d\n', n_increase);

fprintf('\n--- STEP 9 BUILD-PLAN GATES ---\n');
A = r_err_final < 3.0;
fprintf('  (A) Ph4 entry pos err < 3.0 m   : %.3f m   %s\n', r_err_final, pf(A));
B = v_err_final < 0.5;
fprintf('  (B) Ph4 entry vel err < 0.5 m/s : %.3f m/s %s\n', v_err_final, pf(B));
C = max_lcvx_N < 0.1;
fprintf('  (C) LCvx slack < 0.1 N          : %.3e N %s\n', max_lcvx_N, pf(C));
D = mean_iter <= 5;
fprintf('  (D) mean SCP iter <= 5          : %.2f     %s\n', mean_iter, pf(D));
E = ~hit_imax;
fprintf('  (E) no i_max hits               : %s        %s\n', tf(hit_imax), pf(E));
F = sigma_jump < 1.5 && n_increase <= 1;
fprintf('  (F) sigma smooth (jump<1.5, ups<=1) : jump=%.2fs, ups=%d %s\n', ...
        sigma_jump, n_increase, pf(F));
G = max_solve_ms < 50;
fprintf('  (G) max solve time < 50 ms      : %.1f ms %s\n', max_solve_ms, pf(G));

all_pass = A && B && C && D && E && F && G;
fprintf('\nSTEP 9 STRICT PASS : %s\n', pf(all_pass));

%%  PLOTS

figure('Name','Mode B closed loop (Step 7 revert)', ...
       'Position',[100 100 1300 800]);

% 3D trajectory
% 3D trajectory
subplot(2,3,1);
h_traj  = plot3(x_log(1,:), x_log(2,:), x_log(3,:), 'b-', 'LineWidth', 1.4); hold on;
h_start = plot3(x0_si(1), x0_si(2), x0_si(3), 'go', 'MarkerSize', 8, 'MarkerFaceColor','g');
h_target = plot3(target_state(1), target_state(2), target_state(3), ...
                 'r*', 'MarkerSize', 14, 'LineWidth', 1.5);

% Dummy markers — give the legend something to point at for each marker type
h_valid   = plot3(NaN, NaN, NaN, 'k^', 'MarkerSize', 8, 'LineWidth', 1.2);
h_invalid = plot3(NaN, NaN, NaN, 'rx', 'MarkerSize', 8, 'LineWidth', 1.2);

% Real markers
for j = 1:numel(mpc_calls)
    idx = find(t_grid >= mpc_calls(j).t + t_phase3_start, 1);
    if ~isempty(idx)
        col = 'k^'; if ~mpc_calls(j).valid, col = 'rx'; end
        plot3(x_log(1,idx), x_log(2,idx), x_log(3,idx), col, ...
              'MarkerSize', 8, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
grid on; axis equal; xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
legend([h_traj, h_start, h_target, h_valid, h_invalid], ...
       {'trajectory','start','Ph4 entry','MPC valid','MPC invalid'}, ...
       'Location','best');
title('3D trajectory');
% Altitude
subplot(2,3,2);
plot(t_grid - t_phase3_start, x_log(3,:), 'b-', 'LineWidth', 1.4); grid on;
xlabel('t since Ph3 [s]'); ylabel('z [m]');
title('Altitude');

% Top-down view
subplot(2,3,3);
plot(x_log(1,:), x_log(2,:), 'b-', 'LineWidth', 1.4); hold on;
plot(target_state(1), target_state(2), 'r*', 'MarkerSize', 12);
grid on; axis equal; xlabel('X [m]'); ylabel('Y [m]');
title('Top-down view');

% Thrust magnitude
subplot(2,3,4);
plot(t_grid - t_phase3_start, T_cmd_log, 'b-', 'LineWidth', 1.4); hold on;
yline(params.Tmin, 'r:', 'T_{min}'); yline(params.Tmax, 'r:', 'T_{max}');
grid on; xlabel('t since Ph3 [s]'); ylabel('||T_{cmd}|| [N]');
title('Thrust magnitude');

% Speed
subplot(2,3,5);
plot(t_grid - t_phase3_start, vecnorm(x_log(4:6,:)), 'b-', 'LineWidth', 1.4);
grid on; xlabel('t since Ph3 [s]'); ylabel('|v| [m/s]');
title('Speed');

% Per-call diagnostics
subplot(2,3,6);
if ~isempty(mpc_calls)
    yyaxis left;
    bar([mpc_calls.t], [mpc_calls.iter], 'FaceColor',[0.2 0.5 0.9],'FaceAlpha',0.6);
    ylabel('SCP iterations');
    yline(params.scp.i_max,'r:','i_{max}');
    yyaxis right;
    plot([mpc_calls.t], [mpc_calls.lcvx_N], 'rs-', 'LineWidth', 1.4);
    ylabel('LCvx slack [N]');
    yline(0.1, 'r:', 'tol');
end
grid on; xlabel('MPC call t [s]');
title('SCP iter & LCvx per call');

sgtitle(sprintf('Step 9 (descent)  r_{err}=%.2f m  v_{err}=%.2f m/s  fuel=%.3f kg  max-lcvx=%.2e N', ...
    r_err_final, v_err_final, fuel_burned, max_lcvx_N));

figure('Name','Sigma convergence');
plot([mpc_calls.t], [mpc_calls.sigma_sol_s], 'b-o', 'LineWidth', 1.4, 'MarkerSize', 6);
hold on;
yline(params.scp.sigma_desc_low,  'r:', '\sigma_{low}');
yline(params.scp.sigma_desc_high, 'r:', '\sigma_{high}');
grid on; xlabel('MPC call t [s]'); ylabel('\sigma [s]');
title('Sigma convergence across MPC calls');
end




function s = pf(ok)
    if ok, s = '[ PASS ]'; else, s = '[ FAIL ]'; end
end

function s = tf(b)
    if b, s = 'yes'; else, s = 'no'; end
end