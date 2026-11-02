% main_ThreeDimRocket - Main script to solve the Optimal Control Problem
%
% VTVL Rocket (Four-phase Formulation)
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk
%
% Altered by L. White, University of Sheffield May 2026


%% Solve the problem
% main_ThreeDimRocket.m
% Multi-start version for fuel-optimal global optimum search
%
% Strategy: run the solver N_starts times from different initial guesses.
% With 12s per solve, N_starts = 50 takes ~10 minutes.
% If all successful solves converge to the same fuel cost (std dev < 0.05 kg),
% the global optimum is likely found.

clear all; close all; format compact;clc;

% Portable
script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);

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
W_BIG  = 4.65;

phase_colors = [0      0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.4660 0.6740 0.1880;
    0.4940 0.1840 0.5560];
phase_names  = {'Ascent','Hover','Descent','Land'};


%% USER SETTINGS — adjust these
N_starts   = 50;     % number of solver attempts (increase for more confidence)
rng(42);             % fixed random seed - makes results reproducible
                     % change this number to get a different random sequence

%%  PROBLEM SETUP - done once outside the loop
%  (loading problem definition is slow; only do it once)
options.mp = settings_ThreeDimRocketFourPhase;
[problem, guess_original, options.phaseoptions] = ThreeDimRocketFourPhase;

%% Turn off IPOPT console output for all runs except the first
% (avoids hundreds of lines of solver output in your terminal)
options_quiet       = options;
options_quiet.mp.ipopt.print_level = 0;   % 0 = silent, 5 = verbose

%%  STORAGE
best_cost     = inf;        % track minimum fuel used [kg]
best_solution = [];         % store best solution struct
all_fuel      = nan(N_starts, 1);     % fuel used per attempt [kg]
all_feasible  = false(N_starts, 1);   % did the attempt converge?

% Initial mass is fixed — used to compute fuel consumed
m_initial = problem.phases{1}.states.x0(7);   % = ms + mp = 16.88 kg

%%  MULTI-START LOOP
fprintf('Starting multi-start optimisation: %d attempts\n', N_starts);
fprintf('Initial mass: %.3f kg\n\n', m_initial);

for k = 1:N_starts

    fprintf('--- Attempt %d / %d ---  ', k, N_starts);

    % DEBUG - check every phase has constraintErrorTol
    for i = 1:length(problem.phases)
        fprintf('Phase %d: constraintErrorTol exists = %d, settings = %s\n', ...
            i, isfield(problem.phases{i}, 'constraintErrorTol'), ...
            func2str(problem.phases{i}.settings));
    end

    if k == 1
        % First attempt: original guess unchanged - gives a good baseline
        guess_k  = guess_original;
        opts_k   = options;         % verbose output for first solve only
        fprintf('(original guess)\n');
    else
        % All subsequent attempts: randomly perturb the original guess
        % NOTE:  perturb from the ORIGINAL guess, not the last
        % solution. This ensures the search space is fully explored rather
        % than converging around one local minimum.
        guess_k  = perturbGuess(guess_original, problem);
        opts_k   = options_quiet;   % silent output
        fprintf('(perturbed guess)\n');
    end

    try
        [sol_k, ~] = solveMyProblem(problem, guess_k, opts_k);
        fprintf('    Solver returned, n_phases in solution = %d\n', length(sol_k.phaseSol));

        % ---- Compute fuel used [kg] ----
        % Total fuel = initial mass - final mass after phase 4.
        % Mass is state variable x(:,7), continuous across phases via
        % the linkage constraints, so we only need the Phase 4 endpoint.
        n_phases  = length(sol_k.phaseSol);
        sol_last  = sol_k.phaseSol{n_phases};
        m_final   = speval(sol_last, 'X', 7, sol_last.tf);
        fuel_used = m_initial - m_final;

        all_fuel(k)     = fuel_used;
        all_feasible(k) = true;

        fprintf('    Fuel used: %.4f kg\n', fuel_used);

        % Keep the solution with minimum fuel
        if fuel_used < best_cost
            best_cost     = fuel_used;
            best_solution = sol_k;
            fprintf('    *** New best: %.4f kg ***\n', best_cost);
        end

    catch ME
        % Solver failed to converge — skip this attempt
        fprintf('    Solver failed: %s\n', getReport(ME, 'extended'));
    end

end


%%  RESULTS SUMMARY
feasible_fuels = all_fuel(all_feasible);

fprintf('MULTI-START RESULTS SUMMARY\n');
fprintf('-------------------------------------\n');
fprintf('Total attempts:       %d\n',   N_starts);
fprintf('Successful solves:    %d\n',   sum(all_feasible));
fprintf('Failed solves:        %d\n',   sum(~all_feasible));
fprintf('\nFuel consumption across successful solves:\n');
fprintf('  Best  (minimum):    %.4f kg\n', min(feasible_fuels));
fprintf('  Worst (maximum):    %.4f kg\n', max(feasible_fuels));
fprintf('  Mean:               %.4f kg\n', mean(feasible_fuels));
fprintf('  Std deviation:      %.4f kg\n', std(feasible_fuels));
fprintf('\nINTERPRETATION:\n');
fprintf('  Std dev < 0.05 kg  -> global optimum very likely found\n');
fprintf('  Std dev 0.05-0.2   -> increase N_starts for more confidence\n');
fprintf('  Std dev > 0.2      -> multiple local minima exist; increase N_starts\n');
fprintf('--------------------------------------------\n\n');


%%  USE BEST SOLUTION FOR FIGURE GENERATION

solution = best_solution;

% Saves the fuel-optimal result permanently for comparison
save('fuel_optimal_solution.mat', 'best_solution', 'best_cost', 'all_fuel');
fprintf('Fuel-optimal solution saved to fuel_optimal_solution.mat\n');

%% Export reference trajectory for online MPC
% Converts ICLOCS 11-state Euler representation to 13-state quaternion format
% State mapping: [x,y,z,vx,vy,vz,m,theta_deg,phi_deg,p_rad/s,q_rad/s]
%             -> [x,y,z,vx,vy,vz,m,q0,q1,q2,q3,p_rad/s,q_rad/s]
% psi is pinned to 0 throughout (yaw uncontrolled in ICLOCS).

N_export = 1000;   % points per phase
reference.phase_transition_times = zeros(1, length(solution.phaseSol)+1);
reference.t     = [];
reference.x     = [];   % Nx13: quaternion-based
reference.u     = [];   % Nx3:  [T(N), alpha(deg), beta(deg)]
reference.phase = [];

for i = 1:length(solution.phaseSol)
    sol_i = solution.phaseSol{i};
    tt    = linspace(sol_i.t0, sol_i.tf, N_export);

    % --- Extract 11 ICLOCS states ---
    x11 = zeros(N_export, 11);
    for j = 1:11
        x11(:,j) = speval(sol_i, 'X', j, tt)';
    end

    % --- Extract controls [T, alpha_deg, beta_deg] ---
    u3 = zeros(N_export, 3);
    for j = 1:3
        u3(:,j) = speval(sol_i, 'U', j, tt)';
    end

    % --- Convert Euler angles to quaternion (ZYX: yaw=0, pitch=theta, roll=phi) ---
    theta_rad = deg2rad(x11(:,8));   % pitch [rad]
    phi_rad   = deg2rad(x11(:,9));   % roll  [rad]
    psi_rad   = zeros(N_export, 1);  % yaw   [rad] - always 0 in ICLOCS

    % ZYX quaternion: q = qZ * qY * qX
    q0 = cos(psi_rad/2).*cos(theta_rad/2).*cos(phi_rad/2) ...
       + sin(psi_rad/2).*sin(theta_rad/2).*sin(phi_rad/2);
    q1 = cos(psi_rad/2).*cos(theta_rad/2).*sin(phi_rad/2) ...
       - sin(psi_rad/2).*sin(theta_rad/2).*cos(phi_rad/2);
    q2 = cos(psi_rad/2).*sin(theta_rad/2).*cos(phi_rad/2) ...
       + sin(psi_rad/2).*cos(theta_rad/2).*sin(phi_rad/2);
    q3 = sin(psi_rad/2).*cos(theta_rad/2).*cos(phi_rad/2) ...
       - cos(psi_rad/2).*sin(theta_rad/2).*sin(phi_rad/2);

    % Normalise (speval interpolation can drift slightly)
    qnorm = sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    q0 = q0./qnorm;  q1 = q1./qnorm;
    q2 = q2./qnorm;  q3 = q3./qnorm;

    % --- Build 13-state matrix ---
    % [x, y, z, vx, vy, vz, m, q0, q1, q2, q3, p(rad/s), q_rate(rad/s)]
    % Note: ICLOCS states 10,11 (p,q) are in rad/s (SI torque/inertia dynamics)
    x13 = [x11(:,1:7), q0, q1, q2, q3, x11(:,10), x11(:,11)];

    % --- Accumulate ---
    reference.t     = [reference.t;     tt'];
    reference.x     = [reference.x;     x13];
    reference.u     = [reference.u;     u3];
    reference.phase = [reference.phase; i*ones(N_export,1)];

    reference.phase_transition_times(i)   = sol_i.t0;
    reference.phase_transition_times(i+1) = sol_i.tf;
end

reference.fuel_used = best_cost;
reference.params.ms = 8;
reference.params.mp = 8.88;
reference.params.c  = 1488;
reference.params.n_phases = 4;

save('reference_trajectory.mat', 'reference');
fprintf('Reference trajectory exported: %d points, %.4f kg fuel\n', ...
        length(reference.t), reference.fuel_used);

% Quick sanity check
assert(~any(isnan(reference.x(:))), 'NaN in exported states');
assert(~any(isnan(reference.u(:))), 'NaN in exported controls');
assert(all(abs(reference.t - sort(reference.t)) < 1e-9), 'Time not monotonic');
% Check monotonicity within each phase (tight tolerance)
% Phase boundaries can have tiny speval interpolation artefacts — checked separately
for ph = 1:length(solution.phaseSol)
    idx  = reference.phase == ph;
    m_ph = reference.x(idx, 7);
    assert(all(diff(m_ph) <= 1e-3), ...
        sprintf('Mass increasing within phase %d (max rise = %.2e kg)', ...
        ph, max(diff(m_ph))));
end

% Check phase-boundary jumps are not physically meaningful (< 10g is interpolation noise)
boundary_idx = find(diff(reference.phase) ~= 0);
for bi = boundary_idx'
    m_jump = reference.x(bi+1, 7) - reference.x(bi, 7);
    assert(m_jump <= 0.01, ...
        sprintf('Physically large mass jump at phase boundary: +%.4f kg', m_jump));
end
fprintf('Mass monotonicity check passed.\n');
q_norms = sqrt(sum(reference.x(:,8:11).^2, 2));
assert(all(abs(q_norms - 1) < 1e-6), 'Quaternion norm error');
fprintf('Export sanity checks passed.\n');


%% Generate Figures
filepath = 'figures/test_16/';
if ~exist(filepath,'dir'), mkdir(filepath); end
n_phases = length(solution.phaseSol);
%% Figure: 3D trajectory (isometric view)
figure('Name','3D Trajectory', ...
    'Units','inches','Position',[1 1 W_FULL 4.8], ...
    'PaperPositionMode','auto');
hold on;

h_phase = gobjects(n_phases,1);
for i = 1:n_phases
    sol = solution.phaseSol{i};
    tt  = linspace(sol.t0, sol.tf, 500);
    h_phase(i) = plot3(speval(sol,'X',1,tt), speval(sol,'X',2,tt), ...
        speval(sol,'X',3,tt), '-', ...
        'Color', phase_colors(i,:), 'LineWidth', 1.5);
end

xlabel('$X$ [m]'); ylabel('$Y$ [m]'); zlabel('$Z$ [m]');
title('Reference trajectory: 3D view','Interpreter','latex','FontSize',12);
grid on; box on; view(135,15);
xlim([-50 250]);
ylim([-250 50]);
zlim([0 350]);
axis equal;

ax = gca; ax.Position = [0.10, 0.18, 0.85, 0.78];

lgd = legend(h_phase, phase_names, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, lgd.Position(3), lgd.Position(4)];

exportgraphics(gcf, fullfile(filepath,'fig_ref_trajectory_3d.pdf'), ...
    'ContentType','vector');

%% Figure: Position vs time (X, Y, Z stacked)
figure('Name','Position vs Time', ...
    'Units','inches','Position',[1 1 W_FULL 4.6], ...
    'PaperPositionMode','auto');

tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

axis_labels = {'$X$ [m]', '$Y$ [m]', '$Z$ [m]'};
h_phase = gobjects(n_phases,1);

for k = 1:3
    nexttile; hold on;
    for i = 1:n_phases
        sol = solution.phaseSol{i};
        tt  = linspace(sol.t0, sol.tf, 500);
        h = plot(tt, speval(sol,'X',k,tt), '-', ...
            'Color', phase_colors(i,:), 'LineWidth', 1.5);
        if k == 1, h_phase(i) = h; end
    end
    ylabel(axis_labels{k});
    if k == 3
        xlabel('Time [s]');
        ylim([0 400]);
    else
        ylim([-300 300]);
    end
    grid on; box on;
end

title(tl,'Reference position vs time','Interpreter','latex','FontSize',12);

lgd = legend(h_phase, phase_names, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Layout.Tile = 'south';

exportgraphics(gcf, fullfile(filepath,'fig_ref_position.pdf'), ...
    'ContentType','vector');
%% Figure: Velocity vs time (Vx, Vy, Vz stacked)
figure('Name','Velocity vs Time', ...
    'Units','inches','Position',[1 1 W_FULL 4.6], ...
    'PaperPositionMode','auto');

tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

axis_labels = {'$V_x$ [m/s]', '$V_y$ [m/s]', '$V_z$ [m/s]'};
h_phase = gobjects(n_phases,1);

for k = 1:3
    nexttile; hold on;
    for i = 1:n_phases
        sol = solution.phaseSol{i};
        tt  = linspace(sol.t0, sol.tf, 500);
        h = plot(tt, speval(sol,'X',k+3,tt), '-', ...
            'Color', phase_colors(i,:), 'LineWidth', 1.5);
        if k == 1, h_phase(i) = h; end
    end
    ylabel(axis_labels{k});
    if k == 3
        xlabel('Time [s]');
        ylim([-50 50]);
    else
        ylim([-15 15]);
    end
    grid on; box on;
end

title(tl,'Reference velocity vs time','Interpreter','latex','FontSize',12);

lgd = legend(h_phase, phase_names, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Layout.Tile = 'south';

exportgraphics(gcf, fullfile(filepath,'fig_ref_velocity.pdf'), ...
    'ContentType','vector');
%% Figure: Mass vs time
figure('Name','Mass vs Time', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
hold on;

h_phase = gobjects(n_phases,1);
for i = 1:n_phases
    sol = solution.phaseSol{i};
    tt  = linspace(sol.t0, sol.tf, 500);
    h_phase(i) = plot(tt, speval(sol,'X',7,tt), '-', ...
        'Color', phase_colors(i,:), 'LineWidth', 1.5);
end

yline(problem.phases{1}.states.xl(7), '--', 'No fuel line', ...
    'Color',[0.6 0 0], 'Interpreter','latex','FontSize',9, ...
    'LabelVerticalAlignment','bottom','HandleVisibility','off');
yline(problem.phases{1}.states.xu(7), '--', 'Max fuel line', ...
    'Color',[0.6 0 0], 'Interpreter','latex','FontSize',9, ...
    'LabelVerticalAlignment','top','HandleVisibility','off');

xlabel('Time [s]'); ylabel('Mass [kg]');
title('Reference mass vs time','Interpreter','latex','FontSize',12);
ylim([0 18]);
grid on; box on;

ax = gca; ax.Position = [0.10, 0.22, 0.85, 0.68];

lgd = legend(h_phase, phase_names, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, lgd.Position(3), lgd.Position(4)];

exportgraphics(gcf, fullfile(filepath,'fig_ref_mass.pdf'), ...
    'ContentType','vector');
%% Figure: Thrust vs time
figure('Name','Thrust vs Time', ...
    'Units','inches','Position',[1 1 W_FULL 4.0], ...
    'PaperPositionMode','auto');
hold on;

h_phase = gobjects(n_phases,1);
for i = 1:n_phases
    sol = solution.phaseSol{i};
    tt  = linspace(sol.t0, sol.tf, 500);
    h_phase(i) = plot(tt, speval(sol,'U',1,tt), '-', ...
        'Color', phase_colors(i,:), 'LineWidth', 1.5);
end

yline(problem.phases{1}.inputs.uu(1), '--', 'Max thrust', ...
    'Color',[0.6 0 0], 'LineWidth',0.8, ...
    'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left', ...
    'Interpreter','latex', 'FontSize',9, 'HandleVisibility','off');

xlabel('Time [s]'); ylabel('Thrust [N]');
title('Reference thrust vs time','Interpreter','latex','FontSize',12);
grid on; box on;

ax = gca; ax.Position = [0.10, 0.22, 0.85, 0.68];

lgd = legend(h_phase, phase_names, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Units = 'normalized';
lgd.Position = [0.5 - lgd.Position(3)/2, 0.04, lgd.Position(3), lgd.Position(4)];

exportgraphics(gcf, fullfile(filepath,'fig_ref_thrust.pdf'), ...
    'ContentType','vector');

%% Figure: Attitude angles vs time
figure('Name','Attitude vs Time', ...
    'Units','inches','Position',[1 1 W_FULL 4.4], ...
    'PaperPositionMode','auto');

tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

axis_labels = {'$\theta$ [deg]', '$\phi$ [deg]'};
state_idx   = [8, 9];
h_phase = gobjects(n_phases,1);

for k = 1:2
    nexttile; hold on;
    for i = 1:n_phases
        sol = solution.phaseSol{i};
        tt  = linspace(sol.t0, sol.tf, 500);
        h = plot(tt, speval(sol,'X',state_idx(k),tt), '-', ...
            'Color', phase_colors(i,:), 'LineWidth', 1.5);
        if k == 1, h_phase(i) = h; end
    end
    yline(problem.phases{1}.states.xl(state_idx(k)), ':', ...
        'Color',[0.6 0 0], 'HandleVisibility','off');
    yline(problem.phases{1}.states.xu(state_idx(k)), ':', ...
        'Color',[0.6 0 0], 'HandleVisibility','off');
    ylabel(axis_labels{k});
    if k == 2, xlabel('Time [s]'); end
    ylim([-45 45]);
    grid on; box on;
end

title(tl,'Reference attitude vs time','Interpreter','latex','FontSize',12);

lgd = legend(h_phase, phase_names, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Layout.Tile = 'south';

exportgraphics(gcf, fullfile(filepath,'fig_ref_attitude.pdf'), ...
    'ContentType','vector');
%% Figure: Body rates vs time
figure('Name','Body Rates vs Time', ...
    'Units','inches','Position',[1 1 W_FULL 4.4], ...
    'PaperPositionMode','auto');

tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

axis_labels = {'$p$ [deg/s]', '$q$ [deg/s]'};
state_idx   = [10, 11];
h_phase = gobjects(n_phases,1);

for k = 1:2
    nexttile; hold on;
    for i = 1:n_phases
        sol = solution.phaseSol{i};
        tt  = linspace(sol.t0, sol.tf, 500);
        h = plot(tt, speval(sol,'X',state_idx(k),tt), '-', ...
            'Color', phase_colors(i,:), 'LineWidth', 1.5);
        if k == 1, h_phase(i) = h; end
    end
    yline(problem.phases{1}.states.xl(state_idx(k)), ':', ...
        'Color',[0.6 0 0], 'HandleVisibility','off');
    yline(problem.phases{1}.states.xu(state_idx(k)), ':', ...
        'Color',[0.6 0 0], 'HandleVisibility','off');
    ylabel(axis_labels{k});
    if k == 2, xlabel('Time [s]'); end
    grid on; box on;
end

title(tl,'Reference body rates vs time','Interpreter','latex','FontSize',12);

lgd = legend(h_phase, phase_names, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Layout.Tile = 'south';

exportgraphics(gcf, fullfile(filepath,'fig_ref_rates.pdf'), ...
    'ContentType','vector');
%% Figure: Gimbal angles vs time
figure('Name','Gimbal vs Time', ...
    'Units','inches','Position',[1 1 W_FULL 4.4], ...
    'PaperPositionMode','auto');

tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

axis_labels = {'$\alpha$ [deg]', '$\beta$ [deg]'};
input_idx   = [2, 3];
h_phase = gobjects(n_phases,1);

for k = 1:2
    nexttile; hold on;
    for i = 1:n_phases
        sol = solution.phaseSol{i};
        tt  = linspace(sol.t0, sol.tf, 500);
        h = plot(tt, speval(sol,'U',input_idx(k),tt), '-', ...
            'Color', phase_colors(i,:), 'LineWidth', 1.5);
        if k == 1, h_phase(i) = h; end
    end
    yline(problem.phases{1}.inputs.ul(input_idx(k)), ':', ...
        'Color',[0.6 0 0], 'HandleVisibility','off');
    yline(problem.phases{1}.inputs.uu(input_idx(k)), ':', ...
        'Color',[0.6 0 0], 'HandleVisibility','off');
    ylabel(axis_labels{k});
    if k == 2, xlabel('Time [s]'); end
    grid on; box on;
end

title(tl,'Reference gimbal angles vs time','Interpreter','latex','FontSize',12);

lgd = legend(h_phase, phase_names, ...
    'Orientation','horizontal','Box','off','Interpreter','latex');
lgd.Layout.Tile = 'south';

exportgraphics(gcf, fullfile(filepath,'fig_ref_gimbal.pdf'), ...
    'ContentType','vector');

%%  LOCAL FUNCTION: perturbGuess
%  Generates a physically sensible random perturbation of the
%  initial guess. Called for attempts 2 through N_starts.


function guess_out = perturbGuess(guess_in, problem)
% PERTURBGUESS  Create a randomly perturbed initial guess
%
% Perturbation strategy:
%   1. Phase transition times - varied within their allowed bounds
%   2. State trajectory endpoints - varied around the mission waypoints
%   3. Thrust input guesses - randomised within [T_min, 0.8*T_max]
%
% Gimbal angles (alpha, beta) are kept near zero as large angle guesses
% tend to push IPOPT toward infeasible regions.

guess_out = guess_in;   % copy structure, then modify fields
n_phases  = length(problem.phases);



% 1. PERTURB PHASE TRANSITION TIMES
% ----------------------------------------------------------------
%  t = [0,  t2,    t3,      t4,    t5]
%       |   |      |        |      |
%      ph1 end  ph2 end  ph3 end  ph4 end

%  Constraints:
%    t(1) = 0          (fixed)
%    t(2) in [4, 18]   (phase 1 ascent: minimum 4s, maximum 18s)
%    t(3) = t(2) + 5   (hover is exactly 5s — fixed by boundary constraint)
%    t(4) in [t(3)+4, 38]  (phase 3 descent: minimum 4s)
%    t(5) = t(4) + 5   (landing is exactly 5s — fixed by boundary constraint)


t2 = 4 + rand * 14;                     % phase 1 end: 4 to 18 s
t3 = t2 + 5;                             % hover always 5 s
max_t4 = min(38, t3 + 20);              % cap total mission at ~43 s
t4 = t3 + 4 + rand * (max_t4 - t3 - 4); % phase 3: at least 4 s
t5 = t4 + 5;                             % landing always 5 s

guess_out.mp.time = [0, t2, t3, t4, t5];

% Update each phase's individual t0 and tf to match
guess_out.phases{1}.t0 = 0;    guess_out.phases{1}.tf = t2;
guess_out.phases{2}.t0 = t2;   guess_out.phases{2}.tf = t3;
guess_out.phases{3}.t0 = t3;   guess_out.phases{3}.tf = t4;
guess_out.phases{4}.t0 = t4;   guess_out.phases{4}.tf = t5;


% 2. PERTURB STATE TRAJECTORY GUESSES
% ----------------------------------------------------------------
%  guess.phases{i}.states is a [2 x 11] matrix:
%    row 1 = state at phase start
%    row 2 = state at phase end
%
%  Vary the intermediate waypoints around the target positions.
%  Phase 1 start is fixed (launch from rest at origin)
%  Phase 4 end is fixed (landed, zero velocity)
% ----------------------------------------------------------------
m0    = problem.phases{1}.states.x0(7);   % initial mass 16.88 kg
m_min = problem.phases{1}.states.xl(7);   % dry mass 8 kg

% Perturbation scales  [X, Y, Z, Vx, Vy, Vz, m, theta, phi, p, q]
pos_s = 30;   % metres - spatial perturbation
vel_s = 3;    % m/s   - velocity perturbation
ang_s = 8;    % deg   - angle perturbation
rat_s = 2;    % deg/s - angular rate perturbation



% Apex waypoint (originally hover entry at 100, -100, 300)
% Now used as ascent end / descent start
ep_apex = [100 + pos_s*(rand-0.5), ...
    -100 + pos_s*(rand-0.5), ...
    300 + 10*(rand-0.5), ...
    vel_s*(rand-0.5), vel_s*(rand-0.5), vel_s*(rand-0.5), ...
    m0 * (0.65 + 0.25*rand), ...
    ang_s*(rand-0.5), ang_s*(rand-0.5), ...
    rat_s*(rand-0.5), rat_s*(rand-0.5)];

% Pre-landing waypoint (200, -200, 5)
ep_preland = [200 + pos_s*(rand-0.5), ...
    -200 + pos_s*(rand-0.5), ...
    5 + 2*rand, ...
    vel_s*(rand-0.5), vel_s*(rand-0.5), -1*rand, ...
    m0 * (0.40 + 0.25*rand), ...
    ang_s*(rand-0.5), ang_s*(rand-0.5), ...
    rat_s*(rand-0.5), rat_s*(rand-0.5)];

% Hover Disabled
% Phase 1 end / Phase 2 start  — near hover target (100, -100, 300)
ep1 = [100 + pos_s*(rand-0.5), ...
      -100 + pos_s*(rand-0.5), ...
       300 + 10*(rand-0.5), ...          % Z: tighter, near 300 m
       vel_s*(rand-0.5), vel_s*(rand-0.5), vel_s*(rand-0.5), ...
       m0 * (0.65 + 0.25*rand), ...      % mass: 65-90% of initial
       ang_s*(rand-0.5), ang_s*(rand-0.5), ...
       rat_s*(rand-0.5), rat_s*(rand-0.5)];

% Phase 2 end / Phase 3 start  — still near hover (small drift from ep1)
ep2 = ep1 + [5*(rand-0.5), 5*(rand-0.5), 0.3*(rand-0.5), ...
             1*(rand-0.5), 1*(rand-0.5), 1*(rand-0.5), ...
             -0.05*rand, ...             % slight mass decrease during hover
             2*(rand-0.5), 2*(rand-0.5), rat_s*(rand-0.5), rat_s*(rand-0.5)];

% Phase 3 end / Phase 4 start  — near pre-landing target (200, -200, 5)
ep3 = [200 + pos_s*(rand-0.5), ...
      -200 + pos_s*(rand-0.5), ...
         5 + 2*rand, ...                 % Z slightly above 5 m
       vel_s*(rand-0.5), vel_s*(rand-0.5), -1*rand, ... % slight downward velocity
       m0 * (0.40 + 0.25*rand), ...      % mass: 40-65% of initial
       ang_s*(rand-0.5), ang_s*(rand-0.5), ...
       rat_s*(rand-0.5), rat_s*(rand-0.5)];


% original zero-state landing guess
guess_out.phases{1}.states(2,:) = ep1;
guess_out.phases{2}.states(1,:) = ep1;
guess_out.phases{2}.states(2,:) = ep2;
guess_out.phases{3}.states(1,:) = ep2;
guess_out.phases{3}.states(2,:) = ep3;
guess_out.phases{4}.states(1,:) = ep3;

% Clip all mass guesses to physically valid range [m_dry, m_initial]
for i = 1:n_phases
    guess_out.phases{i}.states(:,7) = ...
        max(m_min, min(m0, guess_out.phases{i}.states(:,7)));
end


% 3. PERTURB THRUST INPUT GUESSES
% ----------------------------------------------------------------
%  Thrust is the primary driver of fuel consumption, so varying it
%  widely is important for exploring different trajectory shapes.
%  Angles are kept small to avoid sending IPOPT into bad territory.
% ----------------------------------------------------------------
T_max     = problem.phases{1}.inputs.uu(1);     % 3200 N
T_min_val = problem.phases{1}.inputs.ul(1);     % 32 N  (1% of max)

for i = 1:n_phases
    T_s = T_min_val + rand * (0.75 * T_max - T_min_val);
    T_e = T_min_val + rand * (0.75 * T_max - T_min_val);
    guess_out.phases{i}.inputs(1,1) = T_s;
    guess_out.phases{i}.inputs(2,1) = T_e;

    guess_out.phases{i}.inputs(1,2) = 2*(rand-0.5);   % alpha start
    guess_out.phases{i}.inputs(2,2) = 2*(rand-0.5);   % alpha end
    guess_out.phases{i}.inputs(1,3) = 2*(rand-0.5);   % beta start
    guess_out.phases{i}.inputs(2,3) = 2*(rand-0.5);   % beta end
end

end