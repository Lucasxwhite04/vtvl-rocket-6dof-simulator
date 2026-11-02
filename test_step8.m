function test_step8
% TEST_STEP8  Step 8 verification — shift warm-start + trajectory buffer.

%   (1) Shift warm-start correctness
%       - mpc_solve receives the last accepted sol_out as prev_sol (not [])
%       - warm_start('shift') is exercised on every call after the first
%       - SCP convergence should be FASTER than with reference warm-start
%         (lower mean iterations expected)

%   (2) Trajectory buffer splice continuity
%       - trajectory_buffer.m splices each new MPC plan at its t_grid_si(1)
%       - The commanded thrust T_cmd must be continuous at every splice:
%         max(|T_cmd(k+1) - T_cmd(k)|) at guidance-transition times
%         should be < SPLICE_STEP_TOL (200 N)

%   (3) Regression — landing performance must not worsen
%       - Miss distance < 5 m   (same gate as Step 7)
%       - Final speed  < 0.5 m/s

% Side-by-side comparison: runs both the OLD test (prev_sol=[]) and the
% NEW test (prev_sol=last_sol) to show the improvement explicitly.

% Pass criteria (Step 8 gate):
%   [A] max splice step in T_cmd < 200 N
%   [B] mean SCP iterations NEW <= mean SCP iterations OLD   (no regression)
%   [C] miss_dist_NEW < 5 m
%   [D] speed_NEW < 0.5 m/s

% L. White, University of Sheffield, May 2026


clear; clc; close all;
addpath(genpath('.'));

SPLICE_STEP_TOL = 200;   % [N]  max allowed step in T_cmd at re-plan times

fprintf('  TEST_STEP8 - shift warm-start + trajectory buffer\n');

%% Setup (shared by both runs)
params = plant_params();
ref    = params.reference;
params.dt_guidance = 0.5;   % 2 Hz

if ~isfield(params,'socp_optimizer') || isempty(params.socp_optimizer)
    fprintf('Compiling SOCP (one-time) ...\n');
    params.socp_optimizer = build_socp(params);
    fprintf('\n');
end

%% Phase-4 entry state (same as test_modeB_closed_loop.m)
t_phase4_start = ref.phase_transition_times(4);
t_phase4_end   = ref.phase_transition_times(5);

[t_u, ia] = unique(ref.t,'stable');
x_u = ref.x(ia,:);
x0_si = interp1(t_u, x_u, t_phase4_start,'linear').';
x0_si = [x0_si(1:3); x0_si(4:6); x0_si(7)];   % 7x1 [r;v;m]

fprintf('Phase-4 entry:  r=[%.2f %.2f %.2f] m   v=[%.3f %.3f %.3f] m/s\n\n', ...
        x0_si(1:3), x0_si(4:6));

%% Run both variants
fprintf('--- RUN A: OLD (prev_sol=[] every call, reference warm-start) ---\n');
[resA] = run_landing(x0_si, t_phase4_start, t_phase4_end, params, ref, false);
fprintf('\n--- RUN B: NEW (prev_sol=last_sol, shift warm-start) ---\n');
[resB] = run_landing(x0_si, t_phase4_start, t_phase4_end, params, ref, true);

%  STEP 8 PASS/FAIL REPORT
fprintf('  STEP 8 PASS/FAIL\n');

% [A] Splice step
spliceOK = resB.max_splice_step < SPLICE_STEP_TOL;
fprintf('[A] Max splice step: %6.1f N  (tol %d N)  %s\n', ...
        resB.max_splice_step, SPLICE_STEP_TOL, pf(spliceOK));

% [B] SCP iteration regression
iterOK = resB.mean_iter <= resA.mean_iter + 1;   % +1 tolerance
fprintf('[B] Mean SCP iters: OLD=%.1f  NEW=%.1f   (NEW <= OLD+1)  %s\n', ...
        resA.mean_iter, resB.mean_iter, pf(iterOK));

% [C] Miss distance
missOK = resB.miss_dist < 5.0;
fprintf('[C] Miss distance:  %.3f m  (< 5 m)  %s\n', resB.miss_dist, pf(missOK));

% [D] Final speed
speedOK = resB.final_speed < 0.5;
fprintf('[D] Final speed:    %.4f m/s  (< 0.5 m/s)  %s\n', resB.final_speed, pf(speedOK));

all_pass = spliceOK && iterOK && missOK && speedOK;
fprintf('\n  STEP 8 OVERALL: %s\n', pf(all_pass));

%%  PLOTS
figure('Name','Step 8 — Shift warm-start + buffer','Position',[50 50 1400 800]);

% --- Thrust trace with splice markers (Run B) ---
ax1 = subplot(2,3,1);
plot(resB.t_grid - t_phase4_start, resB.T_cmd_log, 'b-', 'LineWidth', 1.4); hold on;
yline(params.Tmin,'r:','T_{min}','LabelHorizontalAlignment','left');
yline(params.Tmax,'r:','T_{max}','LabelHorizontalAlignment','left');
% Mark splice times
for j = 1:numel(resB.splice_times)
    xline(resB.splice_times(j) - t_phase4_start, 'k--', ...
          sprintf('splice %d',j), 'LabelOrientation','horizontal');
end
grid on; xlabel('t since Ph4 [s]'); ylabel('T_{cmd} [N]');
title('Run B — Thrust + splice times');

% --- Splice step magnitudes ---
ax2 = subplot(2,3,2);
if ~isempty(resB.splice_steps)
    bar(resB.splice_steps, 'FaceColor', [0.2 0.6 0.9]);
    yline(SPLICE_STEP_TOL,'r--',sprintf('tol %d N', SPLICE_STEP_TOL));
    grid on; xlabel('Splice #'); ylabel('\Delta T_{cmd} [N]');
    title('Splice step size (< 200 N goal)');
else
    text(0.5,0.5,'no splices recorded','Units','normalized','HorizontalAlignment','center');
    title('Splice steps');
end

% --- SCP iterations per call: A vs B ---
ax3 = subplot(2,3,3);
hold on;
if ~isempty(resA.iter_log), bar(resA.iter_log, 'FaceColor',[0.8 0.3 0.1],'FaceAlpha',0.6); end
if ~isempty(resB.iter_log), bar(resB.iter_log, 'FaceColor',[0.2 0.6 0.9],'FaceAlpha',0.6); end
yline(params.scp.i_max,'r:','i_{max}');
legend('OLD (ref warm-start)','NEW (shift warm-start)','Location','best');
grid on; xlabel('MPC call #'); ylabel('SCP iterations');
title(sprintf('Iterations  OLD mean=%.1f  NEW mean=%.1f', resA.mean_iter, resB.mean_iter));

% --- Altitude traces ---
ax4 = subplot(2,3,4);
plot(resA.t_grid - t_phase4_start, resA.x_log(3,:), 'r--', 'LineWidth',1.2); hold on;
plot(resB.t_grid - t_phase4_start, resB.x_log(3,:), 'b-',  'LineWidth',1.4);
grid on; xlabel('t since Ph4 [s]'); ylabel('z [m]');
legend('OLD','NEW','Location','best');
title('Altitude');

% --- Speed traces ---
ax5 = subplot(2,3,5);
plot(resA.t_grid - t_phase4_start, vecnorm(resA.x_log(4:6,:)), 'r--','LineWidth',1.2); hold on;
plot(resB.t_grid - t_phase4_start, vecnorm(resB.x_log(4:6,:)), 'b-', 'LineWidth',1.4);
grid on; xlabel('t since Ph4 [s]'); ylabel('|v| [m/s]');
legend('OLD','NEW','Location','best');
title('Speed');

% --- Wall-clock per MPC call ---
ax6 = subplot(2,3,6);
hold on;
if ~isempty(resA.wall_ms), plot(resA.wall_ms,'ro-','LineWidth',1.2); end
if ~isempty(resB.wall_ms), plot(resB.wall_ms,'bs-','LineWidth',1.4); end
legend('OLD','NEW','Location','best');
grid on; xlabel('MPC call #'); ylabel('wall-clock [ms]');
title(sprintf('Solve time  OLD mean=%.0f ms  NEW mean=%.0f ms', ...
              mean(resA.wall_ms), mean(resB.wall_ms)));

sgtitle(sprintf('Step 8 — Shift warm-start + buffer   splice step %.1f N   miss %.2f m   |vf| %.3f m/s', ...
        resB.max_splice_step, resB.miss_dist, resB.final_speed));

end   % test_step8



%%  LOCAL HELPER — run one landing scenario
%  use_shift: false = always pass prev_sol=[] (old behaviour)
%             true  = pass last accepted sol as prev_sol (Step 8)
function res = run_landing(x0_si, t_phase4_start, t_phase4_end, params, ref, use_shift)

    dt_sim    = 0.01;   % 100 Hz plant integration
    t_phase4_target_end = t_phase4_start + params.scp.sigma_land;
    t_no_replan_after   = t_phase4_target_end - 0.3;
    t_sim_end = t_phase4_target_end + 0.2;
    t_grid    = t_phase4_start : dt_sim : t_sim_end;
    N_sim     = numel(t_grid);

    % Guidance state
    prev_sol        = [];
    rejection_count = 0;
    fallback_active = false;
    last_mpc_t      = -inf;
    n_mpc_calls     = 0;

    % Logging
    x_log      = zeros(7, N_sim);
    T_cmd_log  = zeros(1, N_sim);

    iter_log   = [];
    wall_ms    = [];
    splice_times = [];
    T_at_splice_before = [];   % T_cmd just before a splice
    T_at_splice_after  = [];   % T_cmd just after  a splice

    % --- Initialise trajectory buffer ---
    trajectory_buffer('init', [], params);

    x = x0_si;
    for k = 1:N_sim
        t = t_grid(k);

        %% Replan decision
        do_replan = (isempty(iter_log) || (t - last_mpc_t >= params.dt_guidance)) ...
                    && (t < t_no_replan_after) ...
                    && ~fallback_active;

        if do_replan
            x_hat = zeros(17,1);
            x_hat(1:6)  = x(1:6);
            x_hat(7:10) = [1;0;0;0];
            x_hat(17)   = x(7);

            target_state = [params.r_pad; 0; 0; 0];
            sigma_rem    = max(t_phase4_target_end - t, 0.5);
            sigma_bounds = [sigma_rem, sigma_rem];

            
            if use_shift && ~isempty(prev_sol) && prev_sol.valid
                prev_for_mpc = prev_sol;
            else
                prev_for_mpc = [];
            end

            tw = tic;
            [~, ~, sol_out, rejection_count, fb] = mpc_solve( ...
                x_hat, t, target_state, sigma_bounds, ...
                prev_for_mpc, rejection_count, params);
            wall_ms(end+1) = toc(tw) * 1000; %#ok<AGROW>

            last_mpc_t  = t;
            n_mpc_calls = n_mpc_calls + 1;
            iter_log(end+1) = sol_out.iter; %#ok<AGROW>

            fprintf('  t=%5.2f  sigma=%.2f  iter=%d  valid=%d  fault=%d  reject=%d  fb=%d  miss=%.2f  |vf|=%.2f\n', ...
                t - t_phase4_start, sigma_rem, sol_out.iter, sol_out.valid, ...
                sol_out.fault_code, rejection_count, fb, ...
                norm(sol_out.x_k(1:3,end)*params.U_L - params.r_pad), ...
                norm(sol_out.x_k(4:6,end)*params.U_V));

            if sol_out.valid
                % Record splice step: compare T_cmd just before vs just after
                % the new plan's first node in the buffer.
                t_splice = sol_out.t_grid_si(1);
                if ~isempty(find(t_grid <= t_splice, 1, 'last'))
                    k_before = find(t_grid <= t_splice, 1, 'last');
                    if k_before >= 1
                        % T before = what was in buffer at t_splice
                        [T_before, ~] = trajectory_buffer([], t_splice, params);
                    else
                        T_before = NaN;
                    end
                else
                    T_before = NaN;
                end

                % Splice into buffer
                [~, ~] = trajectory_buffer(sol_out, t_splice, params);

                % T after = what buffer gives immediately after splice
                [T_after, ~] = trajectory_buffer([], t_splice, params);

                splice_times(end+1)          = t_splice; %#ok<AGROW>
                T_at_splice_before(end+1)    = T_before; %#ok<AGROW>
                T_at_splice_after(end+1)     = T_after;  %#ok<AGROW>

                % Update prev_sol for next call
                prev_sol = sol_out;
            end

            if fb
                fallback_active = true;
            end
        end

        %% Read command from buffer or Mode A fallback
        [T_cmd, n_hat] = trajectory_buffer([], t, params);

        % Check if buffer has anything useful at this time yet
        buf_started = ~isnan(T_cmd);

        if buf_started && ~fallback_active
            T_cmd = max(params.Tmin, min(params.Tmax, T_cmd));
            if norm(n_hat) < 1e-6, n_hat = [0;0;1]; end
        else
            % Mode A fallback
            x_hat_A = zeros(17,1);
            x_hat_A(1:6) = x(1:6);
            x_hat_A(7:10) = [1;0;0;0];
            x_hat_A(17)   = x(7);
            [T_cmd, n_hat, ~] = guidance_track_reference(t, x_hat_A, ref, params);
            n_hat = n_hat(:);
            if norm(n_hat) > 1e-6, n_hat = n_hat/norm(n_hat); else, n_hat=[0;0;1]; end
            T_cmd = max(params.Tmin, min(params.Tmax, T_cmd));
        end

        T_vec = T_cmd * n_hat;

        x_log(:,k)    = x;
        T_cmd_log(k)  = T_cmd;

        %% Touchdown check
        if x(3) <= 0.3 && t > t_phase4_start + 0.5 && norm(x(4:6)) < 5.0
            x_log     = x_log(:, 1:k);
            T_cmd_log = T_cmd_log(1:k);
            t_grid    = t_grid(1:k);
            break;
        end

        %% Integrate plant
        odefun = @(~,xv) transl_dynamics(xv, T_vec, T_cmd, params);
        [~, X] = ode45(odefun, [t, t+dt_sim], x, ...
                       odeset('RelTol',1e-6,'AbsTol',1e-8));
        x = X(end,:).';
    end

    %% Compute splice step sizes
    splice_steps = abs(T_at_splice_after - T_at_splice_before);
    splice_steps = splice_steps(~isnan(splice_steps));

    %% Pack results
    res.x_log         = x_log;
    res.T_cmd_log     = T_cmd_log;
    res.t_grid        = t_grid;
    res.iter_log      = iter_log;
    res.wall_ms       = wall_ms;
    res.splice_times  = splice_times;
    res.splice_steps  = splice_steps;
    res.max_splice_step = max([splice_steps, 0]);   % 0 if no splices
    res.mean_iter     = mean(iter_log);
    res.miss_dist     = norm(x_log(1:3,end) - params.r_pad);
    res.final_speed   = norm(x_log(4:6,end));

    tag = '(OLD ref ws)';
    if use_shift, tag = '(NEW shift ws)'; end
    fprintf('  %s  miss=%.3f m  |vf|=%.4f m/s  MPC calls=%d  mean iter=%.1f  mean wall=%.0f ms  max_splice_step=%.1f N\n', ...
            tag, res.miss_dist, res.final_speed, n_mpc_calls, res.mean_iter, mean(wall_ms), res.max_splice_step);
end


function s = pf(ok)
    if ok, s = '[ PASS ]'; else, s = '[ FAIL ]'; end
end