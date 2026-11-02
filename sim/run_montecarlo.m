function run_montecarlo
% RUN_MONTECARLO  Monte Carlo dispersion analysis for phase-4 landing.
%
% Runs N_RUNS samples of the Mode B phase-4 closed loop with dispersions
% sampled from the bounds in plant_params(). Produces Fig 4.6, a four-panel
% robustness summary, and saves raw per-run results to montecarlo_results.mat.
%
% Dispersions applied per run (uniform within bounds):
%   - Aerodynamic drag CDxy, CDz   (+/-20%)
%   - Effective exhaust velocity c (+/-2%)
%   - Initial position             (+/-[2.0, 2.0, 1.0] m)
%   - Initial velocity             (+/-[0.3, 0.3, 0.15] m/s)
%   - Initial mass at phase-4 entry (+/-3%)
%
% Only phase 4 is dispersed because:
%   - It is the only phase where Mode B (online SCP-MPC) is active.
%   - Phases 1-3 use Mode A reference tracking which is deterministic.
%
% Outputs:
%   fig4_6_montecarlo.pdf / .png  - 4-panel summary figure
%   montecarlo_results.mat        - raw per-run results
%
% USAGE:
%   run_montecarlo          % uses defaults (N_RUNS = 50)
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
%%  Configuration 
N_RUNS         = 50;     % 50 runs ~ 7 min, 100 ~ 14 min
RNG_SEED_BASE  = 0;      % per-run seed = RNG_SEED_BASE + i


fprintf('  RUN_MONTECARLO - Phase-4 dispersion analysis (N = %d)\n', N_RUNS);
fprintf('------------------------------------------------------------\n\n');

%% Build params + SOCP once (reused across runs) 
fprintf('Compiling SOCP optimiser (one-time) ...\n');
params_nominal = plant_params();
ref = params_nominal.reference;

% Pre-compute reference fuel for fuel-deviation gate
[t_unique, ia] = unique(ref.t, 'stable');
x_unique        = ref.x(ia, :);
t_phase4_start  = ref.phase_transition_times(4);
t_phase4_end    = ref.phase_transition_times(5);
ref_mass_entry  = interp1(t_unique, x_unique(:,7), t_phase4_start, 'linear');
ref_mass_end    = interp1(t_unique, x_unique(:,7), t_phase4_end,   'linear');
ref_fuel_burn   = ref_mass_entry - ref_mass_end;

% Nominal phase-4 entry state
x_at_entry  = interp1(t_unique, x_unique, t_phase4_start, 'linear').';
x0_nominal  = [x_at_entry(1:3); x_at_entry(4:6); x_at_entry(7)];

fprintf('SOCP ready. Reference fuel burn (phase 4 nominal): %.4f kg\n\n', ref_fuel_burn);

%% Pre-allocate result arrays 
miss_dist     = nan(N_RUNS, 1);
final_speed   = nan(N_RUNS, 1);
fuel_burned   = nan(N_RUNS, 1);
fuel_dev_pct  = nan(N_RUNS, 1);
mean_iter_run = nan(N_RUNS, 1);
n_rejected    = nan(N_RUNS, 1);
latched       = false(N_RUNS, 1);
max_lcvx      = nan(N_RUNS, 1);
failed        = false(N_RUNS, 1);
final_pos_xy  = nan(N_RUNS, 2);

%% Outer Monte Carlo loop 
t_mc_start = tic;
print_interval = max(1, floor(N_RUNS/10));
bounds = params_nominal.dispersions.bounds;

for i = 1:N_RUNS
    rng(RNG_SEED_BASE + i);  % reproducible per-run sampling

    % --- Make per-run copy of params and inject dispersions ---
    params = params_nominal;
    params.dispersions.scale.CDxy_scale = 1 + (2*rand-1)*bounds.CDxy_frac;
    params.dispersions.scale.CDz_scale  = 1 + (2*rand-1)*bounds.CDz_frac;
    params.dispersions.scale.c_scale    = 1 + (2*rand-1)*bounds.Isp_frac;

    % Initial state dispersion at phase-4 entry
    dr = [2.0; 2.0; 1.0]  .* (2*rand(3,1) - 1);  % m
    dv = [0.3; 0.3; 0.15] .* (2*rand(3,1) - 1);  % m/s

    x0_si = x0_nominal;
    x0_si(1:3) = x0_si(1:3) + dr;
    x0_si(4:6) = x0_si(4:6) + dv;
    x0_si(7)   = x0_si(7) * (1 + (2*rand-1)*0.03);  % +/-3% initial mass

    % --- Run one phase-4 closed-loop case ---
    try
        result = phase4_closed_loop(params, x0_si, t_phase4_start, ref);
        final_pos_xy(i,:) = result.final_pos(1:2).';
        miss_dist(i)     = result.miss_dist;
        final_speed(i)   = result.final_speed;
        fuel_burned(i)   = result.fuel_burned;
        fuel_dev_pct(i)  = 100*(result.fuel_burned - ref_fuel_burn)/ref_fuel_burn;
        mean_iter_run(i) = result.mean_iter;
        n_rejected(i)    = result.n_rejected;
        latched(i)       = result.latched;
        max_lcvx(i)      = result.max_lcvx;
    catch ME
        failed(i) = true;
        fprintf('  Run %d failed: %s\n', i, ME.message);
    end

    % Progress printout
    if mod(i, print_interval) == 0
        elapsed = toc(t_mc_start);
        eta     = elapsed * (N_RUNS - i) / max(i,1);
        fprintf('  Run %3d/%d  (%3.0f%%)  elapsed %4.0fs  eta %4.0fs\n', ...
                i, N_RUNS, 100*i/N_RUNS, elapsed, eta);
    end
end

t_total = toc(t_mc_start);
fprintf('\nMonte Carlo complete: %d runs in %.1f s (%d failures)\n\n', ...
        N_RUNS, t_total, sum(failed));

%% Summary statistics
ok   = ~failed;
n_ok = sum(ok);

fprintf('--- SUMMARY (N = %d successful runs) ---\n', n_ok);
fprintf('Miss distance:    median = %.3f m, 95th pct = %.3f m, max = %.3f m\n', ...
        median(miss_dist(ok)), prctile(miss_dist(ok), 95), max(miss_dist(ok)));
fprintf('Final speed:      median = %.3f m/s, 95th pct = %.3f m/s, max = %.3f m/s\n', ...
        median(final_speed(ok)), prctile(final_speed(ok), 95), max(final_speed(ok)));
fprintf('Fuel deviation:   median = %+.1f%%, 95th |pct| = %.1f%%, max |pct| = %.1f%%\n', ...
        median(fuel_dev_pct(ok)), prctile(abs(fuel_dev_pct(ok)), 95), max(abs(fuel_dev_pct(ok))));
fprintf('Mean SCP iter:    %.2f\n', mean(mean_iter_run(ok)));
fprintf('Fallback latched: %d/%d (%.1f%%)\n', sum(latched), n_ok, 100*sum(latched)/n_ok);
fprintf('A pass (miss < 0.5 m):           %d/%d (%.1f%%)\n', ...
        sum(miss_dist(ok) < 0.5),  n_ok, 100*sum(miss_dist(ok) < 0.5)/n_ok);
fprintf('B-strict pass (|vf| < 0.1 m/s):  %d/%d (%.1f%%)\n', ...
        sum(final_speed(ok) < 0.1), n_ok, 100*sum(final_speed(ok) < 0.1)/n_ok);
fprintf('B-relaxed pass (|vf| < 0.5 m/s): %d/%d (%.1f%%)\n', ...
        sum(final_speed(ok) < 0.5), n_ok, 100*sum(final_speed(ok) < 0.5)/n_ok);
fprintf('\n');

%%  Fig 4.6: separate dissertation-style figures 

% --- Fig 4.6a — Landing accuracy ---
figure('Name','Monte Carlo - Landing accuracy', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
histogram(miss_dist(ok), 'BinWidth', 0.05, ...
    'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
hold on;
xline(0.5, 'r--', 'gate $= 0.5$ m', ...
    'Interpreter','latex','FontSize',10, 'LineWidth',1.5);
xlabel('Miss distance [m]'); ylabel('Count');
title('Landing accuracy','Interpreter','latex','FontSize',12);
grid on; box on;

ax = gca; ax.Position = [0.10, 0.16, 0.85, 0.76];

exportgraphics(gcf, 'fig4_6a_landing_accuracy.pdf', 'ContentType','vector');

% --- Fig 4.6b — Final velocity ---
figure('Name','Monte Carlo - Final velocity', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
histogram(final_speed(ok), 'BinWidth', 0.05, ...
    'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
hold on;
xline(0.1, 'k--', 'strict 0.1', ...
    'Interpreter','latex','FontSize',9, 'LineWidth',1.2);
xline(0.5, 'r--', 'relaxed 0.5', ...
    'Interpreter','latex','FontSize',10, 'LineWidth',1.5);
xlabel('Final $|v|$ [m/s]'); ylabel('Count');
title('Final velocity','Interpreter','latex','FontSize',12);
grid on; box on;

ax = gca; ax.Position = [0.10, 0.16, 0.85, 0.76];

exportgraphics(gcf, 'fig4_6b_final_velocity.pdf', 'ContentType','vector');

% --- Fig 4.6c — Fuel deviation ---
figure('Name','Monte Carlo - Fuel deviation', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
histogram(fuel_dev_pct(ok), 'BinWidth', 1, ...
    'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
hold on;
xline(0,  'k-',  'ICLOCS', ...
    'Interpreter','latex','FontSize',10, 'LineWidth',1.2);
xline(-5, 'r--', '$\pm 5\%$', ...
    'Interpreter','latex','FontSize',10, 'LineWidth',1.2);
xline(+5, 'r--', '', 'LineWidth',1.2);
xlabel('Fuel deviation from ICLOCS [\%]'); ylabel('Count');
title('Fuel deviation','Interpreter','latex','FontSize',12);
grid on; box on;

ax = gca; ax.Position = [0.10, 0.16, 0.85, 0.76];

exportgraphics(gcf, 'fig4_6c_fuel_deviation.pdf', 'ContentType','vector');

% --- Fig 4.6d — Landing dispersion (top-down) ---
figure('Name','Monte Carlo - Landing dispersion', ...
    'Units','inches','Position',[1 1 W_FULL 4.8], ...
    'PaperPositionMode','auto');

dx = final_pos_xy(ok,1) - params_nominal.r_pad(1);
dy = final_pos_xy(ok,2) - params_nominal.r_pad(2);

scatter(dx, dy, 55, miss_dist(ok), 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',0.6);
hold on;
h_pad = plot(0, 0, 'rp', 'MarkerSize', 18, ...
             'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'LineWidth',1.2);

theta = linspace(0, 2*pi, 200);
h_gate = plot(0.5*cos(theta), 0.5*sin(theta), 'r--', 'LineWidth',1.5);

r95 = prctile(miss_dist(ok), 95);
h_p95 = plot(r95*cos(theta), r95*sin(theta), 'k:', 'LineWidth',1.2);

axis equal; grid on; box on;
xlim_range = max(0.55, 1.2*max(abs([dx; dy])));
xlim([-xlim_range, xlim_range]);
ylim([-xlim_range, xlim_range]);
xlabel('$\Delta X$ from pad [m]'); ylabel('$\Delta Y$ from pad [m]');
title('Landing dispersion (top-down)','Interpreter','latex','FontSize',12);

cb = colorbar;
ylabel(cb, 'Miss distance [m]','Interpreter','latex','FontSize',10);
colormap(gca, 'parula');

ax = gca; ax.Position = [0.12, 0.20, 0.72, 0.72];

lgd = legend([h_pad, h_gate, h_p95], ...
    {'Pad', 'Gate (0.5 m)', sprintf('95th pct (%.2f m)', r95)}, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.03, ...
                lgd.Position(3), lgd.Position(4)];

exportgraphics(gcf, 'fig4_6d_landing_dispersion.pdf', 'ContentType','vector');

fprintf('\nFigures saved:\n');
fprintf('  fig4_6a_landing_accuracy.pdf\n');
fprintf('  fig4_6b_final_velocity.pdf\n');
fprintf('  fig4_6c_fuel_deviation.pdf\n');
fprintf('  fig4_6d_landing_dispersion.pdf\n');

save('montecarlo_results.mat', 'miss_dist', 'final_speed', 'fuel_burned', ...
     'fuel_dev_pct', 'mean_iter_run', 'n_rejected', 'latched', 'max_lcvx', ...
     'failed', 'N_RUNS', 'final_pos_xy');

fprintf('Fig 4.6 exported: fig4_6_montecarlo.{pdf,png}\n');
fprintf('Results saved:    montecarlo_results.mat\n');
fprintf('-------------------------------------------------------\n');

end


function result = phase4_closed_loop(params, x0_si, t_phase4_start, ref)
% PHASE4_CLOSED_LOOP  One phase-4 closed-loop run, returns summary struct.


params.dt_guidance      = 0.5;
dt_sim                  = 0.01;
t_phase4_target_end     = t_phase4_start + params.scp.sigma_land;
t_no_replan_after       = t_phase4_target_end - 0.3;
t_sim_end               = t_phase4_target_end + 0.2;
t_grid                  = t_phase4_start:dt_sim:t_sim_end;
N_sim                   = numel(t_grid);

prev_sol        = [];
rejection_count = 0;
fallback_active = false;
last_mpc_t      = -inf;
buf_t = []; buf_T = []; buf_n = [];

x_log     = zeros(7, N_sim);
mpc_iters = [];
mpc_lcvx  = [];
mpc_valid = [];

x = x0_si;
final_k = N_sim;

for k = 1:N_sim
    t = t_grid(k);

    do_replan = (isempty(buf_t) || (t - last_mpc_t >= params.dt_guidance)) ...
                && (t < t_no_replan_after) && ~fallback_active;

    if do_replan
        x_hat = zeros(17,1);
        x_hat(1:6)  = x(1:6);
        x_hat(7:10) = [1; 0; 0; 0];
        x_hat(17)   = x(7);
        target_state = [params.r_pad; 0; 0; 0];
        sigma_rem    = max(t_phase4_target_end - t, 0.5);
        sigma_bounds = [sigma_rem, sigma_rem];

        [~, ~, sol_out, rejection_count, fb] = mpc_solve( ...
            x_hat, t, target_state, sigma_bounds, ...
            prev_sol, rejection_count, params);
        last_mpc_t = t;

        mpc_iters(end+1) = sol_out.iter;
        mpc_lcvx(end+1)  = sol_out.lcvx_slack_N;
        mpc_valid(end+1) = sol_out.valid;

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
        end
    end

    % Read command
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
    x_log(:, k) = x;

    % Touchdown check
    if x(3) <= 0.01 && t > t_phase4_start + 0.5 && norm(x(4:6)) < 0.3
        final_k = k;
        break;
    end

    % Integrate plant
    odefun = @(~, xv) transl_dynamics(xv, T_vec, Gamma, params);
    [~, X] = ode45(odefun, [t, t+dt_sim], x, ...
                   odeset('RelTol', 1e-6, 'AbsTol', 1e-8));
    x = X(end, :).';
end

x_log = x_log(:, 1:final_k);

% Pack results
result.final_pos   = x_log(1:3, end);
result.miss_dist   = norm(x_log(1:3,end) - params.r_pad);
result.final_speed = norm(x_log(4:6,end));
result.fuel_burned = x0_si(7) - x_log(7,end);
if isempty(mpc_iters)
    result.mean_iter = NaN;
    result.max_lcvx  = NaN;
else
    result.mean_iter = mean(mpc_iters);
    valid_lcvx       = mpc_lcvx(logical(mpc_valid));
    result.max_lcvx  = max(valid_lcvx, [], 'omitnan');
end
result.n_rejected  = rejection_count;
result.latched     = fallback_active;
end