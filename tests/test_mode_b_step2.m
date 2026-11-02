% TEST_MODEB_STEP2  Verify linearise_dynamics.m against nonlinear propagation.
%
% The linearised system is EXACT at the reference iterate. So propagating
% the linear discrete model from x_bar_k with the reference controls and
% sigma must reproduce the result of nonlinearly integrating
% dx/dtau = sigma*f(x, u_aug)  over [tau_k, tau_{k+1}] from the same
% initial state.
%
% If they agree to integration tolerance: linearisation is correct.
% If they disagree: A, B, C, Sigma, or z has a sign/index/units bug.
% L. White, University of Sheffield, May 2026


clear; clc;
addpath(genpath('.'));


fprintf('  TEST_MODEB_STEP2 - linearise_dynamics exactness check\n');
fprintf('--------------------------------------------------------\n\n');

params = plant_params();
ref = params.reference;
N = params.scp.N;

%% Step 1 - build a reference iterate from the landing segment
% Phase 4 starts at t = 33.78 s and ends at t = 36.78 s
t_start = ref.phase_transition_times(4);
t_end   = ref.phase_transition_times(5);
sigma_dim = t_end - t_start;                  % 3.0 s
sigma_nd  = sigma_dim / params.U_T;           % non-dim time dilation

fprintf('Landing segment:  t = [%.3f, %.3f] s   sigma = %.3f s\n', ...
        t_start, t_end, sigma_dim);
fprintf('Non-dim sigma  = %.4f  (= sigma_dim / U_T)\n', sigma_nd);
fprintf('\n');

% De-duplicate phase-boundary repeats (same trick as guidance_track_reference)
[t_unique, ia] = unique(ref.t, 'stable');
x_unique = ref.x(ia, :);
u_unique = ref.u(ia, :);

% Resample to N nodes within the landing segment
t_grid = linspace(t_start, t_end, N).';
x_phase4 = interp1(t_unique, x_unique(:, 1:7), t_grid, 'linear');   % N x 7
u_phase4 = interp1(t_unique, u_unique,         t_grid, 'linear');   % N x 3

% Pre-extract and pre-interpolate the quaternion history
q_all     = x_unique(:, 8:11);                              % 3997 x 4
q_at_grid = interp1(t_unique, q_all, t_grid, 'linear');     % N x 4

% Reconstruct inertial-frame thrust vectors from (T, alpha, beta) and the
% reference quaternion at each node.
T_vec_si = zeros(3, N);
Gamma_si = zeros(1, N);
for k = 1:N
    T_mag = u_phase4(k, 1);
    a_rad = u_phase4(k, 2) * pi/180;
    b_rad = u_phase4(k, 3) * pi/180;
    % Body frame: same sign convention as plant_dynamics Step 4
    n_body = [-sin(a_rad); sin(b_rad); cos(a_rad)*cos(b_rad)];
    q0 = q_at_grid(k,1); q1 = q_at_grid(k,2);
    q2 = q_at_grid(k,3); q3 = q_at_grid(k,4);
    R_bi = [1-2*(q2^2+q3^2),  2*(q1*q2-q0*q3),  2*(q1*q3+q0*q2);
            2*(q1*q2+q0*q3),  1-2*(q1^2+q3^2),  2*(q2*q3-q0*q1);
            2*(q1*q3-q0*q2),  2*(q2*q3+q0*q1),  1-2*(q1^2+q2^2)];
    T_vec_si(:,k) = R_bi * (T_mag * n_body);
    Gamma_si(k)   = T_mag;        % LCvx is tight on the optimum
end

% Convert to non-dim
x_bar     = nondim_state(x_phase4.', params, 'to_nd');     % 7 x N
u_bar     = T_vec_si / params.U_F;                          % 3 x N
Gamma_bar = Gamma_si / params.U_F;                          % 1 x N

fprintf('Reference iterate sampled:  %d nodes\n', N);
fprintf('  x_bar(:,1)   = [%.3f %.3f %.3f | %.3f %.3f %.3f | %.3f]\n', x_bar(:,1));
fprintf('  x_bar(:,end) = [%.3f %.3f %.3f | %.3f %.3f %.3f | %.3f]\n', x_bar(:,end));
fprintf('  u_bar(:,1)   = [%.3f %.3f %.3f]  Gamma_bar(1)   = %.3f\n', ...
        u_bar(:,1), Gamma_bar(1));
fprintf('  u_bar(:,end) = [%.3f %.3f %.3f]  Gamma_bar(end) = %.3f\n', ...
        u_bar(:,end), Gamma_bar(end));
fprintf('\n');

%% Step 2 - run linearise_dynamics
fprintf('Running linearise_dynamics ... ');
tic;
[Abar, Bbar, Cbar, Sigbar, zbar] = ...
    linearise_dynamics(x_bar, u_bar, Gamma_bar, sigma_nd, params);
t_lin = toc;
fprintf('done in %.1f ms\n\n', t_lin * 1000);

%% Step 3 - shape and finite-ness checks
fprintf('Test 3a: output shapes ... ');
assert(numel(Abar) == N-1);
for k = 1:N-1
    assert(isequal(size(Abar{k}),   [7,7]),  'Abar{%d} wrong shape', k);
    assert(isequal(size(Bbar{k}),   [7,4]),  'Bbar{%d} wrong shape', k);
    assert(isequal(size(Cbar{k}),   [7,4]),  'Cbar{%d} wrong shape', k);
    assert(isequal(size(Sigbar{k}), [7,1]),  'Sigbar{%d} wrong shape', k);
    assert(isequal(size(zbar{k}),   [7,1]),  'zbar{%d} wrong shape', k);
end
fprintf('PASS\n');

fprintf('Test 3b: all entries finite ... ');
for k = 1:N-1
    assert(all(isfinite(Abar{k}(:))),   'Abar{%d} has NaN/Inf', k);
    assert(all(isfinite(Bbar{k}(:))),   'Bbar{%d} has NaN/Inf', k);
    assert(all(isfinite(Cbar{k}(:))),   'Cbar{%d} has NaN/Inf', k);
    assert(all(isfinite(Sigbar{k}(:))), 'Sigbar{%d} has NaN/Inf', k);
    assert(all(isfinite(zbar{k}(:))),   'zbar{%d} has NaN/Inf', k);
end
fprintf('PASS\n\n');

%% Step 4 - EXACTNESS CHECK (the main event)
fprintf('Test 4: exactness vs nonlinear integration (the key test)\n');
fprintf('  Per-interval error norm  ||x_lin - x_nl||  for selected k:\n');

dtau = 1/(N-1);
max_err = 0;
for k = 1:N-1
    u_aug_k   = [u_bar(:,k);   Gamma_bar(k)];
    u_aug_kp1 = [u_bar(:,k+1); Gamma_bar(k+1)];
    x_lin = Abar{k} * x_bar(:,k) + ...
            Bbar{k} * u_aug_k + ...
            Cbar{k} * u_aug_kp1 + ...
            Sigbar{k} * sigma_nd + ...
            zbar{k};

    odefun = @(tau, x) sigma_nd * transl_f_nd(tau, x, k, dtau, ...
                          u_bar, Gamma_bar, params);
    [~, X] = ode45(odefun, [(k-1)*dtau, k*dtau], x_bar(:,k), ...
                   odeset('RelTol',1e-10,'AbsTol',1e-12));
    x_nl = X(end, :).';

    err = norm(x_lin - x_nl);
    max_err = max(max_err, err);
    if mod(k, 5) == 0 || k == 1 || k == N-1
        fprintf('    k=%2d:  err = %.3e\n', k, err);
    end
end
fprintf('\n  Max error across all %d intervals: %.3e\n', N-1, max_err);

if max_err < 1e-6
    fprintf('  PASS (max_err < 1e-6)\n');
elseif max_err < 1e-4
    fprintf('  PASS-MARGINAL (max_err < 1e-4) - acceptable\n');
else
    error('FAIL: max_err = %.3e exceeds 1e-4. Bug in linearisation.', max_err);
end
fprintf('\n');

%% Step 5 - 3-DoF-vs-6-DoF mismatch (informational)
fprintf('Test 5: 3-DoF-vs-6-DoF mismatch (informational, not a gate)\n');
mismatch = zeros(N-1, 1);
for k = 1:N-1
    u_aug_k   = [u_bar(:,k);   Gamma_bar(k)];
    u_aug_kp1 = [u_bar(:,k+1); Gamma_bar(k+1)];
    x_lin = Abar{k}*x_bar(:,k) + Bbar{k}*u_aug_k + Cbar{k}*u_aug_kp1 + ...
            Sigbar{k}*sigma_nd + zbar{k};
    mismatch(k) = norm(x_lin - x_bar(:,k+1));
end
fprintf('  Mean mismatch  norm(x_lin - x_bar_{k+1}): %.4e\n', mean(mismatch));
fprintf('  Max  mismatch  norm(x_lin - x_bar_{k+1}): %.4e\n', max(mismatch));
fprintf('  (Reflects how well the 3-DoF planner model matches the 6-DoF\n');
fprintf('   reference. Virtual control nu_k absorbs this in the SCP.)\n\n');

fprintf('--------------------------------------------------------\n');
fprintf('  linearise_dynamics PASSED\n');




%%  LOCAL - nonlinear non-dim dynamics with FOH-interpolated controls
function dx = transl_f_nd(tau, x, k_interval, dtau, u_bar, Gamma_bar, params)
    tau_k = (k_interval - 1) * dtau;
    lam_k   = (tau_k + dtau - tau) / dtau;
    lam_kp1 = (tau - tau_k)        / dtau;
    Tvec = lam_k * u_bar(:,k_interval)   + lam_kp1 * u_bar(:,k_interval+1);
    G    = lam_k * Gamma_bar(k_interval) + lam_kp1 * Gamma_bar(k_interval+1);

    v = x(4:6);
    m = max(x(7), 0.05);
    c_nd = params.c / params.U_V;
    kx = 0.5*params.rho*params.CDxy*params.A_lateral*params.U_V^2/params.U_F;
    kz = 0.5*params.rho*params.CDz *params.A_axial  *params.U_V^2/params.U_F;
    D = [-kx*v(1)*abs(v(1)); -kx*v(2)*abs(v(2)); -kz*v(3)*abs(v(3))];
    g_nd = [0; 0; -1];
    dx = [v; g_nd + Tvec/m + D/m; -G/c_nd];
end