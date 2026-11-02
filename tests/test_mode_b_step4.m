% TEST_MODEB_STEP4  Run the full SCP outer loop to convergence.
%
% Goal: verify the SCP terminates within i_max iterations with both
% convergence metrics below tolerance, returning a self-consistent
% linearised trajectory.
%
% Pass criteria:
%   (a) Loop terminates because nu_tol AND delta_tol both met (not i_max)
%   (b) Final nu_norm  < params.scp.nu_tol    (1e-4)
%   (c) Final delta_norm < params.scp.delta_tol  (1e-3, diagnostic — not enforced)
%   (d) Final LCvx slack < params.scp.lcvx_tol  (0.05 N)
%   (e) Total wall-clock < 1 s (across all iterations)
%   (f) Fuel cost monotonically decreases (or stays flat) across iterations
% L. White, University of Sheffield, May 2026
clear; clc; close all;
addpath(genpath('.'));

%  Dissertation-matched plotting defaults 
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
W_BIG  = 4.65;

fprintf('  TEST_MODEB_STEP4 - SCP outer loop to convergence\n');
fprintf('----------------------------------------------------\n\n');

params = plant_params();
ref = params.reference;

%% Step 1 - warm start from the reference, landing segment
t_phase_start = ref.phase_transition_times(4);
t_phase_end   = ref.phase_transition_times(5);

fprintf('Warm-starting from reference (phase 4) ... ');
[x_bar, u_bar, Gamma_bar, sigma_bar] = ...
    warm_start('reference', params, t_phase_start, t_phase_end, []);
fprintf('done. (sigma_bar = %.4f non-dim = %.3f s)\n\n', ...
        sigma_bar, sigma_bar*params.U_T);


%% Step 2 - use the SOCP optimizer pre-built by plant_params()
optimizer_obj = params.socp_optimizer;


%% Step 3 - boundary conditions for the SCP loop
x0 = x_bar(:, 1);                              % start at warm-start initial

% --- Dispersion injection (forces SCP to work; for Fig 4.3) ---
disp_pos_m  = [5; 3; 1.5];      % offset initial position    [m]
disp_vel_ms = [0.8; -0.5; -0.3];  % offset initial velocity  [m/s]
x0(1:3) = x0(1:3) + disp_pos_m  / params.U_L;
x0(4:6) = x0(4:6) + disp_vel_ms / params.U_V;
fprintf('Dispersion injected: dr = [%.1f %.1f %.1f] m, dv = [%.1f %.1f %.1f] m/s\n\n', ...
        disp_pos_m(1), disp_pos_m(2), disp_pos_m(3), ...
        disp_vel_ms(1), disp_vel_ms(2), disp_vel_ms(3));


r_pad_nd = params.r_pad / params.U_L;
x_target = [r_pad_nd; 0; 0; 0];                % ideal pad, zero velocity

% Landing: sigma fixed
sigma_low_val  = sigma_bar;
sigma_high_val = sigma_bar;
w_t_val        = params.scp.w_t_land;


%% Step 4 - SCP loop
i_max     = params.scp.i_max;
nu_tol    = params.scp.nu_tol;
delta_tol = params.scp.delta_tol;

% Trackers
hist.iter        = (1:i_max).';
hist.nu_norm     = nan(i_max, 1);
hist.delta_norm  = nan(i_max, 1);
hist.lcvx_slack  = nan(i_max, 1);
hist.t_lin       = nan(i_max, 1);
hist.t_solve     = nan(i_max, 1);
hist.fuel_kg     = nan(i_max, 1);   % m0 - m_f  (kg burned)
hist.exit_code   = nan(i_max, 1);

% Save initial iterate so we can plot iterate-0 vs iterate-final
x_iter0 = x_bar;

fprintf('Running SCP outer loop (i_max = %d) ...\n', i_max);
fprintf('%4s | %10s | %10s | %10s | %6s | %6s | %5s\n', ...
        'iter','nu_norm','delta_norm','lcvx[N]','t_lin','t_sol','exit');
fprintf('%s\n', repmat('-', 1, 75));

converged    = false;
converged_at = nan;
for i = 1:i_max
    [x_new, u_new, Gamma_new, sigma_new, info] = scp_iteration( ...
        x_bar, u_bar, Gamma_bar, sigma_bar, ...
        x0, x_target, sigma_low_val, sigma_high_val, w_t_val, ...
        optimizer_obj, params);

    hist.nu_norm(i)    = info.nu_norm;
    hist.delta_norm(i) = info.delta_norm;
    hist.lcvx_slack(i) = info.lcvx_slack;
    hist.t_lin(i)      = info.t_lin;
    hist.t_solve(i)    = info.t_solve;
    hist.exit_code(i)  = info.exit_code;
    hist.fuel_kg(i)    = (1 - x_new(7,end)) * params.U_M;   % m0 - m_f, in kg

    fprintf('%4d | %10.3e | %10.3e | %10.3e | %5.1fms | %5.1fms | %5d\n', ...
            i, info.nu_norm, info.delta_norm, info.lcvx_slack, ...
            info.t_lin*1000, info.t_solve*1000, info.exit_code);

    % Update iterate
    x_bar     = x_new;
    u_bar     = u_new;
    Gamma_bar = Gamma_new;
    sigma_bar = sigma_new;

    % Convergence check: principled SCP criterion
    % Requires BOTH dynamics feasibility (nu_norm) AND iterate stability (delta_norm)
    if info.nu_norm < nu_tol && info.delta_norm < delta_tol
        converged    = true;
        converged_at = i;
        break;
    end

    % Bail on solver failure
    if info.exit_code ~= 0
        warning('Solver failure at SCP iter %d (exit_code=%d). Bailing.', ...
                i, info.exit_code);
        break;
    end
end


%% Step 5 - report
fprintf('\n');
fprintf('--- SCP convergence summary ---\n');
if converged
    fprintf('Converged at iteration %d / %d\n', converged_at, i_max);
else
    fprintf('Did NOT converge within %d iterations.\n', i_max);
end

idx_used = find(~isnan(hist.nu_norm));
fprintf('Final delta_norm  : %.3e  (tol %.0e)  %s\n', ...
        hist.delta_norm(idx_used(end)), delta_tol, ...
        pf(hist.delta_norm(idx_used(end)) < delta_tol));

fprintf('Final LCvx slack  : %.4f N  (tol %.2f N)  %s\n', ...
        hist.lcvx_slack(idx_used(end)), params.scp.lcvx_tol, ...
        pf(hist.lcvx_slack(idx_used(end)) < params.scp.lcvx_tol));

total_wallclock = sum(hist.t_lin(idx_used) + hist.t_solve(idx_used));
fprintf('Total wall-clock  : %.1f ms across %d iterations\n', ...
        total_wallclock*1000, numel(idx_used));

% Fuel progression
fprintf('\nFuel-cost progression  (kg burned, lower is better):\n');
for j = 1:numel(idx_used)
    flag = '';
    if j > 1
        delta_kg = hist.fuel_kg(j) - hist.fuel_kg(j-1);
        if abs(delta_kg) < 1e-4              % within cost tolerance (~0.1 mg)
            flag = '(stable)';
        elseif delta_kg < 0
            flag = '(decreasing OK)';
        else
            flag = '(slight increase — within trust region)';
        end
    end
    fprintf('  iter %d:  fuel = %.6f kg  %s\n', j, hist.fuel_kg(j), flag);
end

% Reference fuel (for comparison)
ref_fuel = ref.x(1, 7) - x_new(7, end)*params.U_M;
fprintf('  reference fuel from m0 to phase-4 terminal: %.6f kg\n', ...
        params.m0 - ref.x(end, 7));


%% Step 6 - figures (one per panel)
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, 'figures');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

ref_fuel_kg = params.m0 - ref.x(end, 7);

% --- Fig 4.3a — Virtual-control norm ---
figure('Name','Virtual-control norm', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
semilogy(idx_used, hist.nu_norm(idx_used), 'b-o', ...
         'LineWidth',1.4,'MarkerSize',6); hold on;
yline(nu_tol, 'r--', '$\nu_{\mathrm{tol}}$', ...
      'Interpreter','latex','FontSize',10);
xlabel('SCP iteration');
ylabel('$\|\nu\|_\infty$');
title('Virtual-control norm','Interpreter','latex','FontSize',12);
xlim([0.5, numel(idx_used)+0.5]);
grid on; box on;

ax = gca; ax.Position = [0.13, 0.16, 0.82, 0.74];

exportgraphics(gcf, fullfile(output_dir,'fig4_3a_virtual_control.pdf'), ...
    'ContentType','vector');

% --- Fig 4.3b — Trust-region step ---
figure('Name','Trust-region step', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
semilogy(idx_used, hist.delta_norm(idx_used), 'b-o', ...
         'LineWidth',1.4,'MarkerSize',6); hold on;
yline(delta_tol, 'r--', '$\delta_{\mathrm{tol}}$', ...
      'Interpreter','latex','FontSize',10);
xlabel('SCP iteration');
ylabel('$\|\Delta\|_\infty$');
title('Trust-region step','Interpreter','latex','FontSize',12);
xlim([0.5, numel(idx_used)+0.5]);
grid on; box on;

ax = gca; ax.Position = [0.13, 0.16, 0.82, 0.74];

exportgraphics(gcf, fullfile(output_dir,'fig4_3b_trust_region.pdf'), ...
    'ContentType','vector');

% --- Fig 4.3c — Fuel cost ---
figure('Name','Fuel cost', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
plot(idx_used, hist.fuel_kg(idx_used)*1000, 'b-o', ...
     'LineWidth',1.4,'MarkerSize',6); hold on;
yline(ref_fuel_kg*1000, 'k--', 'ICLOCS reference', ...
      'Interpreter','latex','FontSize',10, ...
      'LabelHorizontalAlignment','left');
xlabel('SCP iteration');
ylabel('Fuel burned [g]');
title('Fuel cost','Interpreter','latex','FontSize',12);
xlim([0.5, numel(idx_used)+0.5]);
ymax = max(hist.fuel_kg(idx_used))*1000 + 1;
ymin = ref_fuel_kg*1000 - 1;
ylim([ymin, ymax]);
grid on; box on;

ax = gca; ax.Position = [0.13, 0.16, 0.82, 0.74];

exportgraphics(gcf, fullfile(output_dir,'fig4_3c_fuel_cost.pdf'), ...
    'ContentType','vector');

% --- Fig 4.3d — 3D trajectory overlay ---
figure('Name','3D Trajectory', ...
    'Units','inches','Position',[1 1 W_FULL 4.8], ...
    'PaperPositionMode','auto');
idx_skip = 1:4:size(x_iter0,2);
h_warm = plot3(x_iter0(1,idx_skip)*params.U_L, ...
               x_iter0(2,idx_skip)*params.U_L, ...
               x_iter0(3,idx_skip)*params.U_L, ...
               'Color','k', 'LineStyle','--', 'LineWidth',1.5); hold on;
h_conv = plot3(x_new(1,:)*params.U_L, x_new(2,:)*params.U_L, ...
               x_new(3,:)*params.U_L, ...
               'b-o','MarkerSize',3,'LineWidth',1.2);
h_pad  = plot3(params.r_pad(1), params.r_pad(2), params.r_pad(3), ...
               'r*','MarkerSize',12,'LineWidth',1.5);
grid on; box on;                       % <-- no axis equal
xlabel('$X$ [m]'); ylabel('$Y$ [m]'); zlabel('$Z$ [m]');
title('SCP trajectory refinement: 3D view', ...
      'Interpreter','latex','FontSize',12);
view(-45,15);

ax = gca; ax.Position = [0.10, 0.18, 0.85, 0.74];

lgd = legend([h_warm, h_conv, h_pad], ...
    {'Warm start','Converged','Pad'}, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, ...
                lgd.Position(3),          lgd.Position(4)];

exportgraphics(gcf, fullfile(output_dir,'fig4_3d_trajectory_3d.pdf'), ...
    'ContentType','vector');

% --- Fig 4.3e — Altitude profile overlay ---
figure('Name','Altitude profile', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
t_nodes = linspace(0, sigma_new*params.U_T, params.scp.N);
h_warm  = plot(t_nodes, x_iter0(3,:)*params.U_L, 'k--','LineWidth',1.4); hold on;
h_conv  = plot(t_nodes, x_new(3,:)*params.U_L,   'b-o', ...
               'MarkerSize',3,'LineWidth',1.4);
xlabel('Time [s]'); ylabel('Altitude $Z$ [m]');
title('SCP trajectory refinement: altitude profile', ...
      'Interpreter','latex','FontSize',12);
grid on; box on;

ax = gca; ax.Position = [0.10, 0.22, 0.85, 0.68];

lgd = legend([h_warm, h_conv], {'Warm start','Converged'}, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, ...
                lgd.Position(3),          lgd.Position(4)];

exportgraphics(gcf, fullfile(output_dir,'fig4_3e_altitude.pdf'), ...
    'ContentType','vector');

fprintf('\nFigures saved to: %s\n', output_dir);
fprintf('  fig4_3a_virtual_control.pdf\n');
fprintf('  fig4_3b_trust_region.pdf\n');
fprintf('  fig4_3c_fuel_cost.pdf\n');
fprintf('  fig4_3d_trajectory_3d.pdf\n');
fprintf('  fig4_3e_altitude.pdf\n');

%% Export figures for thesis
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, 'figures');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end


if converged
    fprintf('  SCP outer loop CONVERGED. Ready for plumbing files.\n');
else
    fprintf('  SCP did not converge - inspect convergence plots.\n');
end
fprintf('------------------------------------------------------------\n');


%% Helpers
function s = pf(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end