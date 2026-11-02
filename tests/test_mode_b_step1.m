% TEST_MODEB_STEP1  Sanity-check nondim_state.m and transl_dynamics.m.
%
% This is a pure-MATLAB unit test, no Simulink, no SCP yet. Confirms:
%   (1) nondim_state round-trips correctly (to_nd then to_si returns input)
%   (2) Non-dim scales make sense numerically
%   (3) transl_dynamics evaluates without NaN/Inf
%   (4) Trim hover check: at hover thrust, v̇_z ≈ 0
%   (5) Drag sign check: positive vx → negative v̇x contribution
%
% L. White, University of Sheffield, May 2026


clear; clc;
addpath(genpath('.'));


fprintf('  TEST_MODEB_STEP1 — nondim_state + transl_dynamics\n');
fprintf('----------------------------------------------------\n\n');

params = plant_params();

%% Test 1 — round-trip non-dim - SI
fprintf('Test 1: non-dim round trip ... ');
x_si = [100; -50; 200; 5; -3; -10; 14.5];      % typical descent state
x_nd = nondim_state(x_si, params, 'to_nd');
x_back = nondim_state(x_nd, params, 'to_si');
assert(max(abs(x_back - x_si)) < 1e-12, 'Round-trip failed');
fprintf('PASS\n');
fprintf('  x_si    = [%6.1f %6.1f %6.1f | %5.2f %5.2f %5.2f | %5.2f]\n', x_si);
fprintf('  x_nd    = [%6.3f %6.3f %6.3f | %5.3f %5.3f %5.3f | %5.3f]\n', x_nd);
fprintf('\n');

%% Test 2 — non-dim scale sanity
fprintf('Test 2: scale magnitudes ... ');
% Non-dim position should be O(1) — apex z=300 m → 3.0
% Non-dim velocity at terminal ~30 m/s → ~1.0
% Non-dim mass at start = 1.0 (since U_M = m0)
x_apex_nd = nondim_state([0;0;300;0;0;0;16.88], params, 'to_nd');
assert(abs(x_apex_nd(3) - 3.0) < 1e-9, 'Non-dim apex altitude should = 3.0');
m_nd_start = nondim_state([0;0;0;0;0;0;16.88], params, 'to_nd');
assert(abs(m_nd_start(7) - 1.0) < 1e-9, 'Non-dim mass at start should = 1.0');
fprintf('PASS\n');
fprintf('  U_L = %.0f m,  U_T = %.3f s,  U_V = %.2f m/s,  U_M = %.2f kg\n', ...
        params.U_L, params.U_T, params.U_V, params.U_M);
fprintf('  U_F = %.2f N  (expect ≈ 165.5)\n', params.U_F);
fprintf('\n');

%% Test 3 — transl_dynamics evaluates cleanly
fprintf('Test 3: transl_dynamics no NaN/Inf ... ');
x = [0; 0; 200; 0; 0; -10; 14.0];              % mid-descent
T_vec = [0; 0; 140];                            % thrust straight up
Gamma = 140;
xdot = transl_dynamics(x, T_vec, Gamma, params);
assert(all(isfinite(xdot)), 'NaN or Inf in xdot');
fprintf('PASS\n');
fprintf('  xdot = [%.3f %.3f %.3f | %.4f %.4f %.4f | %.5f]\n', xdot);
fprintf('\n');

%% Test 4 — hover trim check
% At hover (v=0), T_vec must balance weight: T/m = g  →  T = m*g
fprintf('Test 4: hover trim (v̇_z ≈ 0 when T = m*g) ... ');
m = 14.0;
T_hover = m * params.g0;                        % expected hover thrust
x = [0; 0; 100; 0; 0; 0; m];
T_vec = [0; 0; T_hover];
xdot = transl_dynamics(x, T_vec, T_hover, params);
% xdot(4:6) should be ≈ 0 since v = 0 (no drag) and T balances g
assert(max(abs(xdot(4:6))) < 1e-9, ...
       'Hover trim failed: max |v̇| = %.3e', max(abs(xdot(4:6))));
fprintf('PASS\n');
fprintf('  T_hover = %.2f N,  v̇ = [%.2e %.2e %.2e]\n', ...
        T_hover, xdot(4), xdot(5), xdot(6));
fprintf('\n');

%% Test 5 — drag sign check
fprintf('Test 5: drag sign (positive v → negative drag accel) ... ');
x = [0; 0; 100; 5; 0; 0; 14.0];                % moving in +x
T_vec = [0; 0; 0];                              % no thrust
Gamma = 0;
xdot = transl_dynamics(x, T_vec, Gamma, params);
% v̇_x should be negative (drag opposes motion) — there's no gravity in x
assert(xdot(4) < 0, 'Drag in +x direction should produce negative v̇_x');
% v̇_z should be -g (gravity only, no vertical drag at v_z = 0)
assert(abs(xdot(6) - params.grav(3)) < 1e-9, 'v̇_z should equal g_z');
fprintf('PASS\n');
fprintf('  v̇_x = %.4f m/s² (drag only)  v̇_z = %.4f m/s² (gravity only)\n', ...
        xdot(4), xdot(6));
fprintf('\n');

%% Done

fprintf('  All 5 tests PASSED. \n');

