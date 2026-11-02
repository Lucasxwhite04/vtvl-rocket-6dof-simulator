%% VERIFY_PLANT.m
% Step 1 verification: plant physics correctness.
%
% Written for the 14-state QUATERNION plant.
%
% State vector (14×1):
% [X; Y; Z; vx; vy; vz; m; q0; q1; q2; q3; p; q_rate; r_rate]
%  1  2  3  4   5   6  7  8   9   10  11  12 13      14
%
% Upright, at-rest initial condition:
% q = [1;0;0;0] (identity rotation — body-z aligned with inertial Z)
% p = q_rate = r_rate = 0
% L. White, University of Sheffield, May 2026

clear; clc;
addpath(fileparts(mfilename('fullpath')));

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

% --- Output folder for thesis figures (absolute path, robust to cwd) ---
script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, 'figures');
if ~exist(output_dir, 'dir')
    [ok, msg] = mkdir(output_dir);
    if ~ok
        error('Could not create %s: %s', output_dir, msg);
    end
end
fprintf('Figures will be saved to: %s\n', output_dir);

% Set a numerical tolerance for floating-point comparisons.
tol = 1e-8;  

params = plant_params();
wind0 = zeros(3, 1);

% --- Convenience: upright, at-rest 14-state initial condition ---
x0_upright = @(z, m) [0; 0; z; 0; 0; 0; m; 1; 0; 0; 0; 0; 0; 0];

%% Test 1 — Hover trim
fprintf('--- Test 1: Hover trim ---\n');
m_hover = params.m0;

% Command a thrust exactly equal to the vehicle's Earth weight.
T_hover = m_hover * params.g0; 
x_hov = x0_upright(100, m_hover);
u_hov = [T_hover; 0; 0];

xd = plant_dynamics(x_hov, u_hov, 0, params, wind0);

% Measure the magnitude (norm) of the resulting rates of change.
pos_vel_err = norm(xd(1:6)); 
quat_rate_err = norm(xd(8:11)); 
ang_acc_err = norm(xd(12:14)); 

fprintf('  xdot(1:6)   = [%.2e, %.2e, %.2e, %.2e, %.2e, %.2e] (pos/vel)\n', xd(1:6));
fprintf('  xdot(7)     = %.4e kg/s (fuel burn -- NOT expected to be 0)\n', xd(7));
fprintf('  xdot(8:11)  = [%.2e, %.2e, %.2e, %.2e] (quat rates -- should be 0)\n', xd(8:11));
fprintf('  xdot(12:14) = [%.2e, %.2e, %.2e] (ang accel -- should be 0)\n', xd(12:14));

if pos_vel_err < tol && quat_rate_err < tol && ang_acc_err < tol
    fprintf('  PASS\n\n');
else
    fprintf('  FAIL\n\n');
end

%% Test 2 — Free-fall acceleration at t=0
fprintf('--- Test 2: Free-fall acceleration at t=0 ---\n');
x_fall = x0_upright(300, params.m0);
u_off = [0; 0; 0]; % Thrust is completely turned off.

xd2 = plant_dynamics(x_fall, u_off, 0, params, wind0);
fprintf('  vz_dot = %.8f (expected %.8f)\n', xd2(6), -params.g0);

if abs(xd2(6) - (-params.g0)) < tol
    fprintf('  PASS\n\n');
else
    fprintf('  FAIL (error = %.2e)\n\n', abs(xd2(6) - (-params.g0)));
end

%% Test 3 — Lateral drag opposes velocity
fprintf('--- Test 3: Lateral drag ---\n');
x_drag = x0_upright(100, params.m0);
x_drag(4) = 10; % Forcibly inject a +10 m/s horizontal velocity.

xd3 = plant_dynamics(x_drag, u_hov, 0, params, wind0);
fprintf('  ax = %.4f m/s^2 (should be < 0 for vx = +10 m/s)\n', xd3(4));
fprintf('  az = %.4f m/s^2 (should be ~0: thrust cancels gravity)\n', xd3(6));

if xd3(4) < 0
    fprintf('  PASS\n\n');
else
    fprintf('  FAIL\n\n');
end

%% Test 4 — Integration: drag SLOWS the fall
fprintf('--- Test 4: Free-fall integration (drag slows fall) ---\n');
x0_fall = x0_upright(300, params.m0);
tspan_dense = linspace(0, 10, 1001); 

%  Use MATLAB's ode45 solver to simulate 5 seconds of falling.
ode_fun = @(t, x) plant_dynamics(x, u_off, t, params, wind0);
[t_out, x_out] = ode45(ode_fun, tspan_dense, x0_fall);

idx_5 = find(abs(t_out - 5.0) < 1e-9, 1);
z_5s = x_out(idx_5, 3); 
z_vacuum = 300 - 0.5 * params.g0 * 5^2; 

% Sets up a simple 1D mathematical model of drag to compare against.
az_1d = @(~, s) [s(2); ...
    -params.g0 - 0.5*params.rho*params.CDz*params.A_axial*s(2)*abs(s(2))/params.m0];
[~, s_1d] = ode45(az_1d, linspace(0, 5, 501), [300; 0]);
z_1d_drag = s_1d(end, 1);

fprintf('  z(t=5 s) 6-DoF sim = %.3f m\n', z_5s);
fprintf('  z_vacuum = %.3f m (lower bound — no drag)\n', z_vacuum);
fprintf('  z_1D_drag_ref = %.3f m (reference — axial drag only)\n', z_1d_drag);

if z_5s > z_vacuum
    fprintf('  PASS: drag correctly slows fall\n\n');
else
    fprintf('  FAIL: z_sim < z_vacuum — drag is acting downward (sign error!)\n\n');
end

%% Test 5 — Gimbal alpha > 0 produces positive pitch acceleration
fprintf('--- Test 5: Gimbal creates pitch acceleration ---\n');
x_att = x0_upright(100, params.m0);

% Command the engine gimbal to pitch (alpha) by +3 degrees.
xd5 = plant_dynamics(x_att, [T_hover; 3*pi/180; 0], 0, params, wind0);
fprintf('  xdot(13) = q_rate_dot = %.6f rad/s^2 (should be > 0)\n', xd5(13));

if xd5(13) > 0
    fprintf('  PASS\n\n');
else
    fprintf('  FAIL\n\n');
end

%% Test 6 — Mass flow rate
fprintf('--- Test 6: Mass flow rate ---\n');

% Verify the mass derivative state (index 7).
xd6 = plant_dynamics(x_hov, u_hov, 0, params, wind0);
mdot_expected = -T_hover / params.c;
fprintf('  mdot = %.6f kg/s (expected %.6f kg/s)\n', xd6(7), mdot_expected);

if abs(xd6(7) - mdot_expected) < tol
    fprintf('  PASS\n\n');
else
    fprintf('  FAIL\n\n');
end

%% Test 7 — Drag symmetry
fprintf('--- Test 7: Drag symmetry ---\n');

% Test drag with the rocket traveling up at 20m/s, then down at 20m/s.
x_up = x0_upright(100, params.m0); x_up(6) = +20;
x_dn = x0_upright(100, params.m0); x_dn(6) = -20;

xd_up = plant_dynamics(x_up, u_hov, 0, params, wind0);
xd_dn = plant_dynamics(x_dn, u_hov, 0, params, wind0);

fprintf('  az(vz=+20 m/s) = %.4f m/s^2 (drag down: should be < 0)\n', xd_up(6));
fprintf('  az(vz=-20 m/s) = %.4f m/s^2 (drag up: should be > 0)\n', xd_dn(6));

if xd_up(6) < 0 && xd_dn(6) > 0
    fprintf('  PASS\n\n');
else
    fprintf('  FAIL\n\n');
end


%%  Separate verification figures (dissertation style)

% --- Shared integrations ---
% Hover (used by Figs 4.2a, 4.2b, 4.2c)
ode_hov = @(t, x) plant_dynamics(x, [x(7)*params.g0; 0; 0], t, params, wind0);
[t_h, xh] = ode45(ode_hov, linspace(0, 5, 501), x_hov);

% Free fall references (used by Fig 4.2d)
ode_vac = @(~, s) [s(2); -params.g0];
[tv, xv] = ode45(ode_vac, linspace(0, 5, 501), [300; 0]);

az_1d = @(~, s) [s(2); ...
    -params.g0 - 0.5*params.rho*params.CDz*params.A_axial*s(2)*abs(s(2))/params.m0];
t_1d = linspace(0, 5, 501);
[~, s_1d_plot] = ode45(az_1d, t_1d, [300; 0]);

mask5     = t_out <= 5.0;
t_6dof_5s = t_out(mask5);
z_6dof_5s = x_out(mask5, 3);

% Drag-model curve (Fig 4.2e)
vz_range = linspace(-50, 50, 300);
Drz = -0.5 * params.rho * params.CDz * params.A_axial * vz_range .* abs(vz_range);

% Gimbal sweep (Fig 4.2f)
alpha_range = linspace(0, 5, 20) * pi/180;
qdot_vals   = zeros(1, numel(alpha_range));
x_att_base  = x0_upright(100, params.m0);
for k = 1:numel(alpha_range)
    xdk = plant_dynamics(x_att_base, [T_hover; alpha_range(k); 0], 0, params, wind0);
    qdot_vals(k) = xdk(13);
end
L_cg2g    = params.ls_nominal;
Iyy_used  = params.Ixy_per_kg * params.m0;
qdot_pred = (T_hover * L_cg2g .* sin(alpha_range)) / Iyy_used;
alpha_deg = alpha_range * 180/pi;

% --- Fig 4.2a — Hover altitude ---
figure('Name','Hover altitude', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
plot(t_h, xh(:,3), 'b-', 'LineWidth', 1.4); hold on;
yline(100, 'r--', 'Reference', 'LabelHorizontalAlignment','right', ...
      'Interpreter','latex','FontSize',9, 'HandleVisibility','off');
xlabel('Time [s]'); ylabel('$z$ [m]');
title('Hover trim: altitude','Interpreter','latex','FontSize',12);
ylim([99.5 100.5]);
grid on; box on;

res_z = max(abs(xh(:,3) - 100));
text(0.05, 0.92, sprintf('$\\max\\,|z-100| = %.2e$ m', res_z), ...
    'Units','normalized','Interpreter','latex','FontSize',10);

ax = gca; ax.Position = [0.10, 0.16, 0.85, 0.76];

exportgraphics(gcf, fullfile(output_dir,'fig4_2a_hover_altitude.pdf'), ...
    'ContentType','vector');

% --- Fig 4.2b — Hover velocity ---
figure('Name','Hover velocity', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
h_vx = plot(t_h, xh(:,4), 'LineWidth', 1.4); hold on;
h_vy = plot(t_h, xh(:,5), 'LineWidth', 1.4);
h_vz = plot(t_h, xh(:,6), 'LineWidth', 1.4);
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Hover trim: velocity','Interpreter','latex','FontSize',12);
ylim([-1 1]);
grid on; box on;

v_max = max(abs([xh(:,4); xh(:,5); xh(:,6)]));
text(0.05, 0.92, sprintf('$\\max\\,|v| = %.2e$ m/s', v_max), ...
    'Units','normalized','Interpreter','latex','FontSize',10);

ax = gca; ax.Position = [0.10, 0.22, 0.85, 0.68];

lgd = legend([h_vx, h_vy, h_vz], {'$v_x$','$v_y$','$v_z$'}, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, lgd.Position(3), lgd.Position(4)];

exportgraphics(gcf, fullfile(output_dir,'fig4_2b_hover_velocity.pdf'), ...
    'ContentType','vector');

% --- Fig 4.2c — Mass depletion ---
figure('Name','Mass depletion', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
h_plant = plot(t_h, xh(:,7), 'b-', 'LineWidth', 1.4); hold on;
mdot_analytical = -T_hover / params.c;
m_pred = xh(1,7) + mdot_analytical * t_h;
h_anal  = plot(t_h, m_pred, 'r--', 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('Mass [kg]');
title('Mass depletion (hover)','Interpreter','latex','FontSize',12);
grid on; box on;

text(0.05, 0.10, sprintf('$\\dot{m} = -T/c = %.4f$ kg/s', mdot_analytical), ...
    'Units','normalized','Interpreter','latex','FontSize',10);

ax = gca; ax.Position = [0.10, 0.22, 0.85, 0.68];

lgd = legend([h_plant, h_anal], {'6-DoF plant','Analytical'}, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, lgd.Position(3), lgd.Position(4)];

exportgraphics(gcf, fullfile(output_dir,'fig4_2c_mass_depletion.pdf'), ...
    'ContentType','vector');

% --- Fig 4.2d - Free fall: drag vs vacuum ---
figure('Name','Free fall', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
h_plant = plot(t_6dof_5s, z_6dof_5s, 'b-', 'LineWidth', 1.4); hold on;
h_vac   = plot(tv, xv(:,1), 'r--', 'LineWidth', 1.0);
xlabel('Time [s]'); ylabel('$z$ [m]');
title('Free fall: drag vs vacuum','Interpreter','latex','FontSize',12);
grid on; box on;

z_6dof_on_1d = interp1(t_6dof_5s, z_6dof_5s, t_1d);
res_drag     = max(abs(z_6dof_on_1d(:) - s_1d_plot(:,1)));
text(0.55, 0.92, sprintf('$|z_{\\mathrm{6DoF}} - z_{\\mathrm{1D}}| < %.1e$ m', res_drag), ...
    'Units','normalized','Interpreter','latex','FontSize',10);

ax = gca; ax.Position = [0.10, 0.22, 0.85, 0.68];

lgd = legend([h_plant, h_vac], {'6-DoF plant','Vacuum'}, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, lgd.Position(3), lgd.Position(4)];

exportgraphics(gcf, fullfile(output_dir,'fig4_2d_free_fall.pdf'), ...
    'ContentType','vector');

% --- Fig 4.2e — Axial drag model ---
figure('Name','Axial drag', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
plot(vz_range, Drz, 'b-', 'LineWidth', 1.4); hold on;
xline(0, '--k', 'HandleVisibility','off');
yline(0, '--k', 'HandleVisibility','off');
xlabel('$v_z$ [m/s]'); ylabel('Axial drag force [N]');
title('Axial drag model','Interpreter','latex','FontSize',12);
grid on; box on;

ax = gca; ax.Position = [0.10, 0.16, 0.85, 0.76];

exportgraphics(gcf, fullfile(output_dir,'fig4_2e_axial_drag.pdf'), ...
    'ContentType','vector');

% --- Fig 4.2f — Pitch acceleration vs gimbal ---
figure('Name','Pitch accel vs gimbal', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
h_plant = plot(alpha_deg, qdot_vals, 'b-',  'LineWidth', 1.4); hold on;
h_anal  = plot(alpha_deg, qdot_pred, 'r--', 'LineWidth', 1.0);
xlabel('Gimbal angle $\alpha$ [deg]');
ylabel('$\dot{q}$ [rad/s$^2$]');
title('Pitch acceleration vs gimbal angle','Interpreter','latex','FontSize',12);
grid on; box on;

ax = gca; ax.Position = [0.10, 0.22, 0.85, 0.68];

lgd = legend([h_plant, h_anal], {'6-DoF plant','Analytical'}, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, lgd.Position(3), lgd.Position(4)];

exportgraphics(gcf, fullfile(output_dir,'fig4_2f_pitch_gimbal.pdf'), ...
    'ContentType','vector');

% Diagnostic: Test 5 used alpha = 3 deg, sweep ends at 5 deg
[~, i3] = min(abs(alpha_deg - 3));
fprintf('\n[Panel 6 diagnostic]\n');
fprintf('  Sweep qdot at alpha = %.2f deg: %.4f rad/s^2\n', alpha_deg(i3), qdot_vals(i3));

fprintf('\nFigures saved to: %s\n', output_dir);
fprintf('  fig4_2a_hover_altitude.pdf\n');
fprintf('  fig4_2b_hover_velocity.pdf\n');
fprintf('  fig4_2c_mass_depletion.pdf\n');
fprintf('  fig4_2d_free_fall.pdf\n');
fprintf('  fig4_2e_axial_drag.pdf\n');
fprintf('  fig4_2f_pitch_gimbal.pdf\n');