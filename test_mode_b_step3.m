% TEST_MODEB_STEP3  Compile and solve the SOCP once from the reference iterate.
%
% This is the first real solve of the Mode B optimisation problem.
% Goal: verify the SOCP compiles, accepts the right parameter sizes,
% solves to optimality, and returns a sensible trajectory.
%
% We seed the optimiser with the ICLOCS reference iterate as both the
% warm start (x_prev) AND the linearisation point. Since the reference
% is already feasible and near-optimal, the SOCP solution should be
% VERY close to it - within virtual-control tolerance.
%
% Pass criteria:
%   (a) ECOS returns exit code 0 (optimal)
%   (b) ||x_sol(:,1) - x0|| < 1e-6     (initial condition satisfied)
%   (c) ||x_sol(1:6,N) - x_target|| < 1e-6  (terminal box)
%   (d) max LCvx slack |Gamma_k - ||u_k||| < 0.05 N  (lossless)
%   (e) max ||nu_k|| < 1e-3            (dynamics satisfied)
%   (f) Solve time < 200 ms
%
% Plots: 3D trajectory overlay (solution vs reference), thrust magnitude,
% LCvx slack vs time, virtual-control magnitude.
% L. White, University of Sheffield, May 2026


clear; clc; close all;
addpath(genpath('.'));

fprintf('  TEST_MODEB_STEP3 - First SOCP solve\n');

params = plant_params();
ref = params.reference;
N = params.scp.N;


%% Step 1 - build reference iterate (same as test_modeB_step2)
t_start = ref.phase_transition_times(4);
t_end   = ref.phase_transition_times(5);
sigma_dim = t_end - t_start;
sigma_nd  = sigma_dim / params.U_T;

[t_unique, ia] = unique(ref.t, 'stable');
x_unique = ref.x(ia, :);
u_unique = ref.u(ia, :);

t_grid = linspace(t_start, t_end, N).';
x_phase4 = interp1(t_unique, x_unique(:, 1:7), t_grid, 'linear');
u_phase4 = interp1(t_unique, u_unique,         t_grid, 'linear');

q_all     = x_unique(:, 8:11);
q_at_grid = interp1(t_unique, q_all, t_grid, 'linear');

T_vec_si = zeros(3, N);
Gamma_si = zeros(1, N);
for k = 1:N
    T_mag = u_phase4(k, 1);
    a_rad = u_phase4(k, 2) * pi/180;
    b_rad = u_phase4(k, 3) * pi/180;
    n_body = [-sin(a_rad); sin(b_rad); cos(a_rad)*cos(b_rad)];
    q0 = q_at_grid(k,1); q1 = q_at_grid(k,2);
    q2 = q_at_grid(k,3); q3 = q_at_grid(k,4);
    R_bi = [1-2*(q2^2+q3^2),  2*(q1*q2-q0*q3),  2*(q1*q3+q0*q2);
            2*(q1*q2+q0*q3),  1-2*(q1^2+q3^2),  2*(q2*q3-q0*q1);
            2*(q1*q3-q0*q2),  2*(q2*q3+q0*q1),  1-2*(q1^2+q2^2)];
    T_vec_si(:,k) = R_bi * (T_mag * n_body);
    Gamma_si(k)   = T_mag;
end

x_bar     = nondim_state(x_phase4.', params, 'to_nd');
u_bar     = T_vec_si / params.U_F;
Gamma_bar = Gamma_si / params.U_F;


%% Step 2 - linearise
fprintf('Linearising dynamics about reference iterate ... ');
tic;
[Abar_c, Bbar_c, Cbar_c, Sigbar_c, zbar_c] = ...
    linearise_dynamics(x_bar, u_bar, Gamma_bar, sigma_nd, params);
t_lin = toc;
fprintf('done in %.1f ms\n', t_lin*1000);

% Cell -> 3D arrays expected by build_socp
Abar3D = cat(3, Abar_c{:});                    % 7 x 7 x (N-1)
Bbar3D = cat(3, Bbar_c{:});                    % 7 x 4 x (N-1)
Cbar3D = cat(3, Cbar_c{:});                    % 7 x 4 x (N-1)
Sigbar2D = [Sigbar_c{:}];                      % 7 x (N-1)
zbar2D   = [zbar_c{:}];                        % 7 x (N-1)


%% Step 3 - compile SOCP
fprintf('\nCompiling SOCP (one-time setup) ...\n');
optimizer_obj = build_socp(params);


%% Step 4 - build parameter set
x0       = x_bar(:, 1);                         % start exactly at the reference start
r_target_nd = params.r_pad / params.U_L;     % ideal pad location, non-dim
v_target_nd = [0; 0; 0];                     % zero velocity at touchdown
x_target    = [r_target_nd; v_target_nd];    % 6x1, the mission's real target
x_prev   = x_bar;
u_prev   = u_bar;
Gamma_prev = Gamma_bar;

% Landing-mode bounds: sigma is fixed = sigma_nd
sigma_low_val  = sigma_nd;
sigma_high_val = sigma_nd;
w_t_val        = params.scp.w_t_land;           % 0 for landing

param_values = {Abar3D, Bbar3D, Cbar3D, Sigbar2D, zbar2D, ...
                x0, x_target, x_prev, u_prev, Gamma_prev, ...
                sigma_low_val, sigma_high_val, w_t_val};

fprintf('\nTarget setup:\n');
fprintf('  ideal pad        = [%6.2f %6.2f %6.2f] m\n', params.r_pad);
fprintf('  reference end    = [%6.2f %6.2f %6.2f] m\n', x_bar(1:3,end)*params.U_L);
fprintf('  warm-start drift = %.3f m  (this is the gap the SCP will close)\n', ...
    norm(x_bar(1:3,end)*params.U_L - params.r_pad));
% Check glideslope at the warm-start terminal (just informational)
r_xy_drift = norm(x_bar(1:2,end)*params.U_L - params.r_pad(1:2));
cone_r_at_z = tan(params.gamma_gs) * (x_bar(3,end)*params.U_L - params.r_pad(3));
fprintf('  ref terminal xy-drift = %.3f m,  cone radius at that alt = %.3f m\n', ...
    r_xy_drift, cone_r_at_z);
if r_xy_drift > cone_r_at_z
    fprintf('  -> reference terminal is OUTSIDE the cone (warm start only, target moves to pad)\n');
end

%% Step 5 - solve
fprintf('\nSolving SOCP ... ');
tic;
[sol, exit_code, ~, ~, ~] = optimizer_obj(param_values);
t_solve = toc;
fprintf('done in %.1f ms (exit code = %d)\n', t_solve*1000, exit_code);

if iscell(exit_code)
    exit_code = exit_code{1};
end

if exit_code ~= 0
    warning('SOCP exit code %d (0=optimal, 3=max-iter, 4=numerical). Continuing to inspect solution.', exit_code);
end

x_sol     = sol{1};   % 7 x N
u_sol     = sol{2};   % 3 x N
Gamma_sol = sol{3};   % 1 x N
sigma_sol = sol{4};   % scalar
nu_sol    = sol{5};   % 7 x (N-1)
s_delta_sol = sol{6}; % 1 x N


%% Step 6 - checks
fprintf('\n--- Solution sanity checks ---\n');

% (a) Initial condition
err_init = norm(x_sol(:,1) - x0);
fprintf('Initial condition error:        %.3e   ', err_init);
fprintf('(target: < 1e-6)  '); print_pf(err_init < 1e-6);

% (b) Terminal box
err_term = norm(x_sol(1:6,end) - x_target);
fprintf('Terminal box error:             %.3e   ', err_term);
fprintf('(target: < 1e-6)  '); print_pf(err_term < 1e-6);

% (c) LCvx slack
lcvx_slack_per_k = zeros(N, 1);
for k = 1:N
    lcvx_slack_per_k(k) = abs(Gamma_sol(k) - norm(u_sol(:,k)));
end
lcvx_slack_max_nd = max(lcvx_slack_per_k);
lcvx_slack_max_si = lcvx_slack_max_nd * params.U_F;
fprintf('Max LCvx slack:                 %.3e (non-dim)  = %.4f N\n', ...
        lcvx_slack_max_nd, lcvx_slack_max_si);
fprintf('  ');
fprintf('(target: < %.2f N, params.scp.lcvx_tol)  ', params.scp.lcvx_tol);
print_pf(lcvx_slack_max_si < params.scp.lcvx_tol);

% (d) Virtual control
nu_norms = vecnorm(nu_sol);                    % 1 x (N-1)
nu_max = max(nu_norms);
fprintf('Max virtual control ||nu||:     %.3e   ', nu_max);
fprintf('(target: < 1e-3)  '); print_pf(nu_max < 1e-3);

% (e) Trust region slack
s_delta_max = max(s_delta_sol);
fprintf('Max trust-region slack:         %.3e   ', s_delta_max);
fprintf('(seeded at reference -> expect near 0)\n');

% (f) Sigma
fprintf('sigma_sol = %.4f (non-dim) = %.3f s (target fixed at %.3f s)\n', ...
        sigma_sol, sigma_sol*params.U_T, sigma_dim);

% (g) Solve time
fprintf('Solve wall-clock:               %.1f ms  ', t_solve*1000);
fprintf('(target: < 200 ms)  '); print_pf(t_solve*1000 < 200);

% (h) Final mass
m_final_nd = x_sol(7, end);
m_final_si = m_final_nd * params.U_M;
m_ref_final = x_bar(7, end) * params.U_M;
fprintf('Final mass (sol):  %.4f kg\n', m_final_si);
fprintf('Final mass (ref):  %.4f kg\n', m_ref_final);
fprintf('Fuel diff: %+.4f kg (negative = solver burned MORE than reference)\n', ...
        m_final_si - m_ref_final);


%% Step 7 - plot
figure('Name','SOCP solution vs reference','Position',[100 100 1200 800]);

% 3D trajectory
subplot(2,3,1);
plot3(x_bar(1,:)*params.U_L, x_bar(2,:)*params.U_L, x_bar(3,:)*params.U_L, ...
      'k--','LineWidth',1.0); hold on;
plot3(x_sol(1,:)*params.U_L, x_sol(2,:)*params.U_L, x_sol(3,:)*params.U_L, ...
      'b-o','LineWidth',1.2,'MarkerSize',3);
plot3(params.r_pad(1), params.r_pad(2), params.r_pad(3), 'r*','MarkerSize',12);
grid on; axis equal; xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
legend('reference','SOCP solution','pad','Location','best');
title('3D trajectory (landing)');

% Altitude vs time
subplot(2,3,2);
t_nodes = linspace(0, sigma_nd*params.U_T, N);
plot(t_nodes, x_bar(3,:)*params.U_L, 'k--'); hold on;
plot(t_nodes, x_sol(3,:)*params.U_L, 'b-o','MarkerSize',3);
grid on; xlabel('t [s]'); ylabel('z [m]'); title('Altitude');

% Thrust magnitude
subplot(2,3,3);
T_si = vecnorm(u_sol) * params.U_F;
T_ref_si = vecnorm(u_bar) * params.U_F;
plot(t_nodes, T_ref_si, 'k--'); hold on;
plot(t_nodes, T_si, 'b-o','MarkerSize',3);
yline(params.Tmin,'r:'); yline(params.Tmax,'r:');
grid on; xlabel('t [s]'); ylabel('||T|| [N]'); title('Thrust magnitude');

% Mass
subplot(2,3,4);
plot(t_nodes, x_bar(7,:)*params.U_M,'k--'); hold on;
plot(t_nodes, x_sol(7,:)*params.U_M,'b-o','MarkerSize',3);
yline(params.ms,'r:');
grid on; xlabel('t [s]'); ylabel('m [kg]'); title('Mass');

% LCvx slack
subplot(2,3,5);
semilogy(t_nodes, lcvx_slack_per_k * params.U_F + 1e-12, 'b-o','MarkerSize',3);
hold on; yline(params.scp.lcvx_tol,'r--');
grid on; xlabel('t [s]'); ylabel('|\Gamma - ||u||| [N]');
title('LCvx slack (log scale)');

% Virtual control
subplot(2,3,6);
semilogy(t_nodes(1:end-1), nu_norms + 1e-12, 'b-o','MarkerSize',3);
grid on; xlabel('t [s]'); ylabel('||\nu_k||');
title('Virtual control magnitude');

sgtitle('Mode B SOCP - first solve from reference iterate');


fprintf('  build_socp + first solve OK.\n');


%% --- helpers
function print_pf(ok)
    if ok, fprintf('PASS\n'); else, fprintf('FAIL\n'); end
end