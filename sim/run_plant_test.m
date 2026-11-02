% RUN_PLANT_TEST  Quick open-loop plant test — gravity freefall for 1 second.
% Verifies plant_dynamics works before building the full closed loop.
% L. White, University of Sheffield, May 2026


clear; clc;
addpath(genpath('..'));

params = plant_params();

% Initial state: on the pad, upright, full mass
x0 = [0;0;0;        % position  (m)
    0;0;0;        % velocity  (m/s)
    16.88;        % mass      (kg)
    1;0;0;0;      % quaternion (upright)
    0;0;0];       % angular rates (rad/s)

% Zero control and wind
u    = zeros(3,1);
wind = zeros(3,1);

% Integrate with ode45 — same solver the build plan specifies
odefun = @(t,x) plant_dynamics(x, u, t, params, wind);
[t_out, x_out] = ode45(odefun, [0 1], x0);

% Check results
fprintf('--- Plant Test Results ---\n');
fprintf('Final z:  %.3f m   (expected ≈ -4.9 m)\n',  x_out(end,3));
fprintf('Final vz: %.3f m/s (expected ≈ -9.8 m/s)\n', x_out(end,6));
fprintf('Final mass: %.4f kg (expected = 16.88 kg, no thrust)\n', x_out(end,7));
fprintf('Quaternion norm: %.8f (expected = 1.0)\n', ...
    norm(x_out(end,8:11)));

% Simple plot
figure;
plot(t_out, x_out(:,3));
xlabel('Time (s)'); ylabel('Altitude z (m)');
title('Freefall test — gravity only');
grid on;