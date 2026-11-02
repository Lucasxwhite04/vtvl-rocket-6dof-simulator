% TEST_MODEB_STEP5  End-to-end Mode B test - single mpc_solve call.
%
% Constructs a realistic x_hat at the start of phase 4 (sampled from
% the locked ICLOCS reference), calls mpc_solve, and verifies it returns
% a usable command + valid sol_out struct.
%
% This is the final algorithmic gate before Simulink integration.
%
% Pass criteria:
%   (a) sol_out.exit_code == 0 (every SCP iteration succeeded)
%   (b) sol_out.valid == true (validity check passed)
%   (c) rejection_count_new == 0 (counter reset on accept)
%   (d) fallback_flag == false
%   (e) T_cmd in [Tmin, Tmax]
%   (f) norm(n_hat_cmd) ≈ 1
%   (g) Total wall-clock < 1 s
%   (h) Final LCvx slack < lcvx_tol
% L. White, University of Sheffield, May 2026


clear; clc; close all;
addpath(genpath('.'));

fprintf('  TEST_MODEB_STEP5 - end-to-end mpc_solve\n');
params = plant_params();
ref = params.reference;

% Compile the SOCP once and stash on params
fprintf('Compiling SOCP (one-time setup) ...\n');
params.socp_optimizer = build_socp(params);
fprintf('\n');


%% Step 1 - construct a realistic x_hat at start of phase 4
% Pull the truth state from the reference at the phase-4 entry time,
% then pad to 17 elements with zero biases (we are not testing the EKF
% here, only the guidance pipeline).
t_now = ref.phase_transition_times(4);                  % start of landing
[t_unique, ia] = unique(ref.t, 'stable');
x_unique = ref.x(ia, :);

% Truth state (13-state plant) at t_now
x_truth = interp1(t_unique, x_unique, t_now, 'linear').';

% Build 17-state x_hat:  [r; v; q; b_g; b_a; m]
x_hat = zeros(17, 1);
x_hat(1:3)   = x_truth(1:3);      % position
x_hat(4:6)   = x_truth(4:6);      % velocity
x_hat(7:10)  = x_truth(8:11);     % quaternion (cols 8-11 of plant state)
x_hat(11:13) = 0;                 % gyro bias (zero for noiseless test)
x_hat(14:16) = 0;                 % accel bias
x_hat(17)    = x_truth(7);        % mass (col 7 of plant state)

fprintf('Scenario:\n');
fprintf('  t_now    = %.3f s\n', t_now);
fprintf('  r_hat    = [%6.2f %6.2f %6.2f] m\n', x_hat(1:3));
fprintf('  v_hat    = [%6.3f %6.3f %6.3f] m/s\n', x_hat(4:6));
fprintf('  m_hat    = %.3f kg\n', x_hat(17));


%% Step 2 - target and time bounds for landing
target_state    = [params.r_pad; 0; 0; 0];                % pad, zero velocity
sigma_bounds_s  = [params.scp.sigma_land, params.scp.sigma_land];  % fixed
prev_sol        = [];                                      % first call
rejection_count = 0;

fprintf('  target   = [%6.2f %6.2f %6.2f] m (pad)\n', target_state(1:3));
fprintf('  sigma_bd = [%.2f %.2f] s\n\n', sigma_bounds_s);


%% Step 3 - call mpc_solve
fprintf('Calling mpc_solve ...\n');
tic;
[T_cmd, n_hat_cmd, sol_out, rejection_count_new, fallback_flag] = ...
    mpc_solve(x_hat, t_now, target_state, sigma_bounds_s, ...
              prev_sol, rejection_count, params);
t_total = toc;
fprintf('done in %.1f ms (wall-clock incl test overhead)\n\n', t_total*1000);


%% Step 4 - report
fprintf('--- mpc_solve outputs ---\n');
fprintf('T_cmd            : %.2f N        (bounds [%.0f, %.0f])\n', ...
        T_cmd, params.Tmin, params.Tmax);
fprintf('||n_hat_cmd||    : %.6f\n', norm(n_hat_cmd));
fprintf('n_hat_cmd        : [%.4f %.4f %.4f]\n', n_hat_cmd);
fprintf('\n');
fprintf('--- sol_out diagnostics ---\n');
fprintf('phase_id         : %d\n',          sol_out.phase_id);
fprintf('SCP iter used    : %d / %d\n',     sol_out.iter, params.scp.i_max);
fprintf('converged early  : %s\n',          tf2str(sol_out.converged));
fprintf('exit_code        : %d (0=optimal)\n', sol_out.exit_code);
fprintf('LCvx slack       : %.4e N  (tol %.2f N)\n', ...
        sol_out.lcvx_slack_N, params.scp.lcvx_tol);
fprintf('nu_norm          : %.4e\n',        sol_out.nu_norm);
fprintf('sigma            : %.4f non-dim = %.3f s\n', ...
        sol_out.sigma, sol_out.sigma*params.U_T);
fprintf('Total t_solve    : %.1f ms\n',     sol_out.t_solve_total*1000);
fprintf('valid            : %s\n',          tf2str(sol_out.valid));
fprintf('fault_code       : %d (0=accept)\n', sol_out.fault_code);
fprintf('rejection_count  : %d -> %d\n', rejection_count, rejection_count_new);
fprintf('fallback_flag    : %s\n',          tf2str(fallback_flag));
fprintf('\n');


%% Step 5 - per-iteration breakdown
fprintf('--- SCP iteration history ---\n');
fprintf('%4s | %10s | %10s | %10s | %10s | %6s | %6s\n', ...
        'iter','nu_norm','delta_norm','lcvx[N]','fuel[kg]','t_lin','t_sol');
fprintf('%s\n', repmat('-', 1, 76));
for k = 1:numel(sol_out.iter_history)
    h = sol_out.iter_history(k);
    fprintf('%4d | %10.3e | %10.3e | %10.3e | %10.6f | %5.1fms | %5.1fms\n', ...
        k, h.nu_norm, h.delta_norm, h.lcvx_slack, h.fuel_kg, ...
        h.t_lin*1000, h.t_solve*1000);
end
fprintf('\n');


%% Step 6 - pass/fail summary
fprintf('--- Pass/fail ---\n');
print_pf('exit_code == 0',              sol_out.exit_code == 0);
print_pf('valid == true',               sol_out.valid);
print_pf('rejection_count_new == 0',    rejection_count_new == 0);
print_pf('fallback_flag == false',      ~fallback_flag);
print_pf('Tmin <= T_cmd <= Tmax',       T_cmd >= params.Tmin && T_cmd <= params.Tmax);
print_pf('||n_hat_cmd|| approx 1',      abs(norm(n_hat_cmd) - 1) < 1e-6);
print_pf('total wall-clock < 1 s',      sol_out.t_solve_total < 1.0);
print_pf('LCvx slack < lcvx_tol',       sol_out.lcvx_slack_N < params.scp.lcvx_tol);


%% Step 7 - plot the trajectory mpc_solve produced
figure('Name','mpc_solve output','Position',[100 100 1200 420]);

subplot(1,3,1);
x_si = sol_out.x_k * params.U_L;                  % approx for plotting
% (Each row scales differently, so do it properly:)
x_dim = zeros(7, params.scp.N);
for k = 1:params.scp.N
    x_dim(:,k) = nondim_state(sol_out.x_k(:,k), params, 'to_si');
end
plot3(x_dim(1,:), x_dim(2,:), x_dim(3,:), 'b-o','MarkerSize',3,'LineWidth',1.2); hold on;
plot3(params.r_pad(1), params.r_pad(2), params.r_pad(3), 'r*','MarkerSize',12);
grid on; axis equal; xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Planned trajectory');

subplot(1,3,2);
plot(sol_out.t_grid_si, x_dim(3,:), 'b-o','MarkerSize',3); 
grid on; xlabel('t [s]'); ylabel('z [m]'); title('Altitude');

subplot(1,3,3);
plot(sol_out.t_grid_si, sol_out.T_k_si, 'b-o','MarkerSize',3); hold on;
yline(params.Tmin,'r:'); yline(params.Tmax,'r:');
grid on; xlabel('t [s]'); ylabel('||T|| [N]'); title('Thrust magnitude');

sgtitle(sprintf('mpc_solve output - first call at t=%.2fs', t_now), ...
        'Interpreter', 'none');


fprintf('  Mode B core algorithm: COMPLETE.\n');

%% --- helpers
function s = tf2str(b), if b, s = 'true'; else, s = 'false'; end, end
function print_pf(label, ok)
    if ok, fprintf('  [ PASS ]  %s\n', label);
    else,  fprintf('  [ FAIL ]  %s\n', label); end
end