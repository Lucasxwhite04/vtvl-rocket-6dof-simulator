function test_modeB_closed_loop
% TEST_MODEB_CLOSED_LOOP  Closed-loop simulation of the Mode B landing phase.
%
% Runs phase 4 (the landing manoeuvre) in closed loop without Simulink.
%   - Initial state taken from the phase-4 entry node of the reference trajectory
%   - 3-DoF translational dynamics integrated by ode45
%   - Mode B MPC replans every params.dt_guidance seconds
%   - Trajectory buffered between replans
%   - Perfect thrust-direction tracking assumed (no attitude loop, no EKF)
%
% Pass criteria (Step 7):
%   (A) Landing position error < 0.5 m
%   (B) Final velocity         < 0.1 m/s
%   (C) LCvx slack             < 0.1 N throughout
%   (D) Mean SCP iterations    <= 5
%   (E) No cycle hits i_max=8
%   (F) Fuel within 5% of ICLOCS reference
%
% Diagnostic output per MPC call:
%   t, sigma, iter, valid, fault_code, exit_code, lcvx_slack_N,
%   pos_err, vel_err, rejection_count, fallback_flag
%
% L. White, University of Sheffield, May 2026
 
clear; clc; close all;
addpath(genpath('.'));

% Dissertation-matched plotting defaults 
set(0, 'DefaultTextInterpreter',          'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter',        'latex');
set(0, 'DefaultAxesFontSize',     11);
set(0, 'DefaultTextFontSize',     11);
set(0, 'DefaultLegendFontSize',   10);
set(0, 'DefaultAxesLineWidth',    0.75);
set(0, 'DefaultLineLineWidth',    1.2);
set(0, 'DefaultLegendBox',        'off');
set(0, 'DefaultAxesBox',          'on');
set(0, 'DefaultFigureColor',      'white');

W_FULL = 6.20;


% ---- Output folder for thesis figures (absolute path) ----
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, 'figures');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
fprintf('Figures will be saved to: %s\n', output_dir);

fprintf('  TEST_MODEB_CLOSED_LOOP - Step 7 \n');

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

%% Initial state: phase-4 entry, sampled from reference
t_phase4_start = ref.phase_transition_times(4);
t_phase4_end   = ref.phase_transition_times(5);

[t_unique, ia] = unique(ref.t, 'stable');
x_unique = ref.x(ia, :);

x_at_entry = interp1(t_unique, x_unique, t_phase4_start, 'linear').';
x0_si = [x_at_entry(1:3); x_at_entry(4:6); x_at_entry(7)];

fprintf('Phase-4 entry conditions:\n');
fprintf('  t_phase4_start = %.3f s\n', t_phase4_start);
fprintf('  r  = [%6.2f %6.2f %6.2f] m\n', x0_si(1:3));
fprintf('  v  = [%6.3f %6.3f %6.3f] m/s\n', x0_si(4:6));
fprintf('  m  = %.3f kg\n', x0_si(7));
fprintf('  target pad = [%6.2f %6.2f %6.2f] m\n\n', params.r_pad);

%% Simulation parameters
dt_sim = 0.01;

t_phase4_target_end = t_phase4_start + params.scp.sigma_land;
fprintf('Landing deadline (fixed): t = %.3f s\n\n', t_phase4_target_end);

t_no_replan_after = t_phase4_target_end - 0.3;
t_sim_end  = t_phase4_target_end + 0.2;
t_grid     = t_phase4_start:dt_sim:t_sim_end;
N_sim      = numel(t_grid);

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
                   'sigma_rem',{}, 'pos_err',{}, 'vel_err',{}, ...
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

        target_state = [params.r_pad; 0; 0; 0];
        sigma_rem    = max(t_phase4_target_end - t, 0.5);
        sigma_bounds = [sigma_rem, sigma_rem];

        t_call = tic;
        [~, ~, sol_out, rejection_count, fb] = mpc_solve( ...
            x_hat, t, target_state, sigma_bounds, ...
            prev_sol, rejection_count, params);
        t_call_ms = toc(t_call) * 1000;

        last_mpc_t  = t;
        n_mpc_calls = n_mpc_calls + 1;

        % Compute terminal errors the SCP is promising
        pos_err_m   = norm(sol_out.x_k(1:3,end)*params.U_L - params.r_pad);
        vel_err_mps = norm(sol_out.x_k(4:6,end)*params.U_V);

        % Pull final-iter exit code from history if available
        exit_code = NaN;
        if isfield(sol_out, 'iter_history') && ~isempty(sol_out.iter_history)
            exit_code = sol_out.iter_history(end).exit_code;
        end

        if sol_out.iter >= params.scp.i_max
            hit_imax = true;
        end

        mpc_calls(end+1) = struct( ...
            't',          t - t_phase4_start, ...
            'valid',      sol_out.valid, ...
            'fault',      sol_out.fault_code, ...
            'exit_code',  exit_code, ...
            'iter',       sol_out.iter, ...
            'wall_ms',    t_call_ms, ...
            'lcvx_N',     sol_out.lcvx_slack_N, ...
            'sigma_rem',  sigma_rem, ...
            'pos_err',    pos_err_m, ...
            'vel_err',    vel_err_mps, ...
            'reject',     rejection_count, ...
            'fb',         fb); %#ok<AGROW>

        fprintf('  %5.2f  %5.2f  %4d  %5d %5d %4d  %7.2e %7.3f  %6.3f  %3d  %2d\n', ...
            t - t_phase4_start, sigma_rem, sol_out.iter, sol_out.valid, ...
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
                t - t_phase4_start);
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

    %% Touchdown check
    if x(3) <= 0.01 && t > t_phase4_start + 0.5 && norm(x(4:6)) < 0.3
        fprintf('  TOUCHDOWN at t=%.2f s\n', t - t_phase4_start);
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

%%  RESULT SUMMARY
final_pos   = x_log(1:3, end);
final_vel   = x_log(4:6, end);
final_mass  = x_log(7, end);
miss_dist   = norm(final_pos - params.r_pad);
final_speed = norm(final_vel);
fuel_burned = x0_si(7) - final_mass;

% Reference fuel for fuel-comparison gate (F)
ref_mass_at_entry = interp1(t_unique, x_unique(:,7), t_phase4_start, 'linear');
ref_mass_at_end   = interp1(t_unique, x_unique(:,7), t_phase4_end,   'linear');
ref_fuel_burn     = ref_mass_at_entry - ref_mass_at_end;

if ~isempty(mpc_calls)
    valid_calls = [mpc_calls.valid];
    if any(valid_calls)
        max_lcvx_N = max([mpc_calls(valid_calls).lcvx_N]);
    else
        max_lcvx_N = NaN;
    end
    mean_iter  = mean([mpc_calls.iter]);
    n_valid    = sum([mpc_calls.valid]);
else
    max_lcvx_N = NaN; mean_iter = NaN; n_valid = 0;
end

fprintf('\n--- LANDING SUMMARY ---\n');
fprintf('Final r          : [%7.3f %7.3f %7.3f] m\n', final_pos);
fprintf('Pad              : [%7.3f %7.3f %7.3f] m\n', params.r_pad);
fprintf('Miss distance    : %.4f m\n',      miss_dist);
fprintf('Final |v|        : %.4f m/s\n',    final_speed);
fprintf('Fuel burned      : %.4f kg  (ICLOCS reference: %.4f kg)\n', ...
        fuel_burned, ref_fuel_burn);
fprintf('Fuel deviation   : %.1f%%\n', 100*(fuel_burned - ref_fuel_burn)/ref_fuel_burn);
fprintf('Total MPC calls  : %d  (%d valid, %d invalid)\n', ...
        n_mpc_calls, n_valid, n_mpc_calls - n_valid);
fprintf('Mean SCP iters   : %.2f\n', mean_iter);
fprintf('Hit i_max=%d     : %s\n', params.scp.i_max, tf(hit_imax));
fprintf('Max LCvx slack   : %.3e N\n', max_lcvx_N);
fprintf('Fallback latched : %s\n', tf(fallback_active));
fprintf('Wall-clock total : %.1f s\n', t_wallclock);

fprintf('\n--- STEP 7 BUILD-PLAN GATES ---\n');
A = miss_dist  < 0.5;        fprintf('  (A) miss < 0.5 m         : %.4f m   %s\n', miss_dist, pf(A));
B = final_speed < 0.1;       fprintf('  (B) |vf| < 0.1 m/s       : %.4f m/s %s\n', final_speed, pf(B));
C = max_lcvx_N < 0.1;        fprintf('  (C) LCvx slack < 0.1 N   : %.3e N %s\n', max_lcvx_N, pf(C));
D = mean_iter <= 5;          fprintf('  (D) mean iter <= 5       : %.2f     %s\n', mean_iter, pf(D));
E = ~hit_imax;               fprintf('  (E) no i_max hits        : %s        %s\n', tf(hit_imax), pf(E));
F = abs(fuel_burned-ref_fuel_burn)/ref_fuel_burn < 0.05;
fprintf('  (F) fuel within 5%%       : %+.1f%%   %s\n', ...
        100*(fuel_burned-ref_fuel_burn)/ref_fuel_burn, pf(F));
fprintf('\n');

% Relaxed gate from project session summary
B_relaxed = final_speed < 0.5;
fprintf('  (B-relaxed) |vf| < 0.5 m/s : %.4f m/s %s   (session-summary gate)\n', ...
        final_speed, pf(B_relaxed));

all_strict = A && B && C && D && E && F;
all_relax  = A && B_relaxed && C && D && E;
fprintf('\nSTRICT BUILD-PLAN PASS : %s\n', pf(all_strict));
fprintf('RELAXED   PASS         : %s\n', pf(all_relax));

%%  Separate figures (dissertation style)

% --- Fig 4.4a — 3D trajectory ---
figure('Name','Phase 4 - 3D trajectory', ...
    'Units','inches','Position',[1 1 W_FULL 4.8], ...
    'PaperPositionMode','auto');
h_traj  = plot3(x_log(1,:), x_log(2,:), x_log(3,:), 'b-', 'LineWidth', 1.4); hold on;
h_start = plot3(x0_si(1), x0_si(2), x0_si(3), 'go', ...
                'MarkerSize', 8, 'MarkerFaceColor','g');
h_pad   = plot3(params.r_pad(1), params.r_pad(2), params.r_pad(3), ...
                'r*', 'MarkerSize', 14, 'LineWidth', 1.5);
h_valid   = plot3(NaN, NaN, NaN, 'k^', 'MarkerSize', 8, 'LineWidth', 1.2);
h_invalid = plot3(NaN, NaN, NaN, 'rx', 'MarkerSize', 8, 'LineWidth', 1.2);

for j = 1:numel(mpc_calls)
    idx = find(t_grid >= mpc_calls(j).t + t_phase4_start, 1);
    if ~isempty(idx)
        col = 'k^'; if ~mpc_calls(j).valid, col = 'rx'; end
        plot3(x_log(1,idx), x_log(2,idx), x_log(3,idx), col, ...
              'MarkerSize', 8, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end

grid on; box on;
xlabel('$X$ [m]'); ylabel('$Y$ [m]'); zlabel('$Z$ [m]');
title('Phase 4 landing: 3D trajectory','Interpreter','latex','FontSize',12);
daspect([1 1 6]);   % preserve original axis stretch
view(40, 20);

ax = gca; ax.Position = [0.10, 0.18, 0.85, 0.74];

lgd = legend([h_traj, h_start, h_pad, h_valid, h_invalid], ...
    {'Trajectory','Start','Pad','MPC valid','MPC invalid'}, ...
    'Orientation','horizontal','Box','off','Interpreter','latex', ...
    'NumColumns',5);
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, ...
                lgd.Position(3), lgd.Position(4)];

exportgraphics(gcf, fullfile(output_dir,'fig4_4a_trajectory_3d.pdf'), ...
    'ContentType','vector');

% --- Fig 4.4b — Altitude ---
figure('Name','Phase 4 - Altitude', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
plot(t_grid - t_phase4_start, x_log(3,:), 'b-', 'LineWidth', 1.4);
grid on; box on;
xlabel('Time since Phase 4 entry [s]'); ylabel('Altitude $z$ [m]');
title('Phase 4 landing: altitude','Interpreter','latex','FontSize',12);

ax = gca; ax.Position = [0.10, 0.16, 0.85, 0.76];

exportgraphics(gcf, fullfile(output_dir,'fig4_4b_altitude.pdf'), ...
    'ContentType','vector');

% --- Fig 4.4c — Top-down view ---
figure('Name','Phase 4 - Top-down', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
plot(x_log(1,:), x_log(2,:), 'b-', 'LineWidth', 1.4); hold on;
plot(params.r_pad(1), params.r_pad(2), 'r*', ...
     'MarkerSize', 12, 'LineWidth', 1.5);
grid on; box on; axis equal;
xlabel('$X$ [m]'); ylabel('$Y$ [m]');
title('Phase 4 landing: top-down view','Interpreter','latex','FontSize',12);

exportgraphics(gcf, fullfile(output_dir,'fig4_4c_topdown.pdf'), ...
    'ContentType','vector');

% --- Fig 4.4d — Thrust magnitude ---
figure('Name','Phase 4 - Thrust', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
plot(t_grid - t_phase4_start, T_cmd_log, 'b-', 'LineWidth', 1.4); hold on;
yline(params.Tmin, 'r:', '$T_{\min}$', ...
      'Interpreter','latex','FontSize',9, ...
      'LabelHorizontalAlignment','left');
yline(params.Tmax, 'r:', '$T_{\max}$', ...
      'Interpreter','latex','FontSize',9, ...
      'LabelHorizontalAlignment','left');
grid on; box on;
xlabel('Time since Phase 4 entry [s]'); ylabel('$\|T_{\mathrm{cmd}}\|$ [N]');
title('Phase 4 landing: thrust magnitude','Interpreter','latex','FontSize',12);

ax = gca; ax.Position = [0.10, 0.16, 0.85, 0.76];

exportgraphics(gcf, fullfile(output_dir,'fig4_4d_thrust.pdf'), ...
    'ContentType','vector');

% --- Fig 4.4e — Speed ---
figure('Name','Phase 4 - Speed', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
plot(t_grid - t_phase4_start, vecnorm(x_log(4:6,:)), 'b-', 'LineWidth', 1.4);
grid on; box on;
xlabel('Time since Phase 4 entry [s]'); ylabel('$|\mathbf{v}|$ [m/s]');
title('Phase 4 landing: speed','Interpreter','latex','FontSize',12);

ax = gca; ax.Position = [0.10, 0.16, 0.85, 0.76];

exportgraphics(gcf, fullfile(output_dir,'fig4_4e_speed.pdf'), ...
    'ContentType','vector');

% --- Fig 4.4f — SCP iterations + LCvx slack ---
if ~isempty(mpc_calls)
    figure('Name','Phase 4 - SCP diagnostics', ...
        'Units','inches','Position',[1 1 W_FULL 4.0], ...
        'PaperPositionMode','auto');

    yyaxis left;
    h_iter = bar([mpc_calls.t], [mpc_calls.iter], ...
        'FaceColor',[0.2 0.5 0.9],'FaceAlpha',0.6,'EdgeColor','none');
    ylabel('SCP iterations');
    yline(params.scp.i_max,'r:','$i_{\max}$', ...
          'Interpreter','latex','FontSize',9);
    ax_l = gca; ax_l.YColor = [0.2 0.5 0.9];

    yyaxis right;
    h_lcvx = plot([mpc_calls.t], [mpc_calls.lcvx_N], 'rs-', 'LineWidth', 1.2);
    ylabel('LCvx slack [N]');
    yline(0.1, 'r:', 'tol', 'Interpreter','latex','FontSize',9);
    ax_r = gca; ax_r.YColor = [0.8 0 0];

    grid on; box on;
    xlabel('MPC call $t$ [s]');
    title('SCP iterations and LCvx slack per call', ...
          'Interpreter','latex','FontSize',12);

    ax = gca; ax.Position = [0.10, 0.22, 0.80, 0.70];

    lgd = legend([h_iter, h_lcvx], {'SCP iterations','LCvx slack'}, ...
        'Orientation','horizontal','Box','off','Interpreter','latex');
    lgd.Units = 'normalized';
    lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, ...
                    lgd.Position(3), lgd.Position(4)];

    exportgraphics(gcf, fullfile(output_dir,'fig4_4f_scp_diagnostics.pdf'), ...
        'ContentType','vector');
end

% --- Fig 4.5 — Solver wall-clock histogram ---
if ~isempty(mpc_calls)
    figure('Name','SCP wall-clock distribution', ...
        'Units','inches','Position',[1 1 W_FULL 4.0], ...
        'PaperPositionMode','auto');

    wall_ms_all = [mpc_calls.wall_ms];
    budget_ms   = params.dt_guidance * 1000;

    histogram(wall_ms_all, 'BinWidth', 25, ...
        'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
    hold on;
    xline(budget_ms, 'r--', sprintf('Budget = %d ms', budget_ms), ...
        'LineWidth', 1.5, 'Interpreter','latex','FontSize',10, ...
        'LabelVerticalAlignment','top');
    xline(median(wall_ms_all), 'k-', ...
        sprintf('Median = %.0f ms', median(wall_ms_all)), ...
        'LineWidth', 1.2, 'Interpreter','latex','FontSize',10);

    xlabel('SCP solve wall-clock [ms]'); ylabel('Count');
    title('SCP-MPC solve time distribution (Phase 4)', ...
          'Interpreter','latex','FontSize',12);
    grid on; box on;

    pct95 = prctile(wall_ms_all, 95);
    pct99 = prctile(wall_ms_all, 99);
    n_overrun = sum(wall_ms_all > budget_ms);
    txt = sprintf(['$n = %d$ calls\n95th pct = %.0f ms\n' ...
                   '99th pct = %.0f ms\nOverruns = %d (%.1f\\%%)'], ...
        numel(wall_ms_all), pct95, pct99, n_overrun, ...
        100*n_overrun/numel(wall_ms_all));
    text(0.97, 0.97, txt, 'Units','normalized', ...
        'HorizontalAlignment','right','VerticalAlignment','top', ...
        'Interpreter','latex','FontSize',10, ...
        'BackgroundColor','w','EdgeColor',[0.7 0.7 0.7],'Margin',6);

    ax = gca; ax.Position = [0.10, 0.16, 0.85, 0.76];

    exportgraphics(gcf, fullfile(output_dir,'fig4_5_solver_timing.pdf'), ...
        'ContentType','vector');
end

fprintf('\nFigures saved to: %s\n', output_dir);
fprintf('  fig4_4a_trajectory_3d.pdf\n');
fprintf('  fig4_4b_altitude.pdf\n');
fprintf('  fig4_4c_topdown.pdf\n');
fprintf('  fig4_4d_thrust.pdf\n');
fprintf('  fig4_4e_speed.pdf\n');
fprintf('  fig4_4f_scp_diagnostics.pdf\n');
fprintf('  fig4_5_solver_timing.pdf\n');

end

function s = pf(ok)
    if ok, s = '[ PASS ]'; else, s = '[ FAIL ]'; end
end

function s = tf(b)
    if b, s = 'yes'; else, s = 'no'; end
end