%  DIAGNOSE_TRACKING.M
%  Run AFTER a simulation to find why the vehicle doesn't follow
%  the reference trajectory.
%
%  Prerequisites: params and reference must already be in workspace,
%  and the simulation must have produced logged signals (logsout or
%  individual workspace variables).


%%  SECTION 1: Reference integrity check 
fprintf('\n--- SECTION 1: Reference Integrity ---\n');

if ~exist('reference','var')
    error('reference not in workspace. Run plant_params() or load reference_trajectory.mat manually.');
end
if ~exist('params','var')
    error('params not in workspace. Run: params = plant_params()');
end

fprintf('reference.t   : %d points,  t=[%.2f ... %.2f] s\n', ...
    numel(reference.t), reference.t(1), reference.t(end));
fprintf('reference.x   : %dx%d\n', size(reference.x));
fprintf('reference.u   : %dx%d\n', size(reference.u));
fprintf('reference.phase: unique phases = %s\n', ...
    mat2str(unique(reference.phase)'));

% Check for duplicate timestamps (phase boundaries)
[~, ia] = unique(reference.t, 'stable');
n_dups  = numel(reference.t) - numel(ia);
fprintf('Duplicate timestamps at phase boundaries: %d  (expected 3)\n', n_dups);

% Check thrust column
T_ref_all = reference.u(:,1);
fprintf('T_ref range: [%.1f, %.1f] N  (Tmin=%.0f, Tmax=%.0f)\n', ...
    min(T_ref_all), max(T_ref_all), params.Tmin, params.Tmax);
if min(T_ref_all) < params.Tmin*0.9
    fprintf('  *** WARNING: T_ref dips below Tmin — reference may be infeasible ***\n');
end

% Check gimbal angles
alpha_all = reference.u(:,2);  beta_all  = reference.u(:,3);
fprintf('alpha_ref range: [%.2f, %.2f] deg  (limit ±%.1f deg)\n', ...
    min(alpha_all), max(alpha_all), params.alpha_max*180/pi);
fprintf('beta_ref  range: [%.2f, %.2f] deg  (limit ±%.1f deg)\n', ...
    min(beta_all),  max(beta_all),  params.beta_max *180/pi);

%%  SECTION 2: Reference trajectory overview plot 
fprintf('\n--- SECTION 2: Reference Overview ---\n');

[t_c, ia] = unique(reference.t, 'stable');
x_c       = reference.x(ia,:);
u_c       = reference.u(ia,:);
ph_c      = reference.phase(ia);

figure('Name','Reference Overview','NumberTitle','off');

subplot(2,3,1);
plot3(x_c(:,1), x_c(:,2), x_c(:,3));
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('3-D Reference Path'); grid on; axis equal;

subplot(2,3,2);
plot(t_c, x_c(:,3));
xlabel('t (s)'); ylabel('Z (m)'); title('Altitude'); grid on;

subplot(2,3,3);
plot(t_c, sqrt(sum(x_c(:,4:6).^2, 2)));
xlabel('t (s)'); ylabel('Speed (m/s)'); title('Reference Speed'); grid on;

subplot(2,3,4);
plot(t_c, u_c(:,1));
xlabel('t (s)'); ylabel('T (N)'); title('Reference Thrust'); grid on;
yline(params.Tmin,'r--','Tmin'); yline(params.Tmax,'r--','Tmax');

subplot(2,3,5);
plot(t_c, u_c(:,2), t_c, u_c(:,3));
legend('\alpha_{ref}','\beta_{ref}'); xlabel('t (s)'); ylabel('deg');
title('Reference Gimbal Angles'); grid on;
yline( params.alpha_max*180/pi,'r--'); yline(-params.alpha_max*180/pi,'r--');

subplot(2,3,6);
plot(t_c, double(ph_c));
xlabel('t (s)'); ylabel('Phase'); title('Phase Index'); grid on;
yticks(1:4); ylim([0.5 4.5]);

%% SECTION 3: Probe guidance outputs at key time points 
fprintf('\n--- SECTION 3: Guidance Spot-Check ---\n');

t_probe = [2, 8, 20, 40, 55];    % seconds - adjust
fprintf('%-6s  %-8s  %-30s  %-8s  %-8s\n', ...
    'Time','T_cmd(N)','n_hat_cmd','pos_err','vel_err');

for k = 1:numel(t_probe)
    tp = t_probe(k);
    if tp > t_c(end), continue; end

    % Build a perfect MEKF state from the reference (no errors)
    r_p = interp1(t_c, x_c(:,1:3), tp, 'linear')';
    v_p = interp1(t_c, x_c(:,4:6), tp, 'linear')';
    m_p = interp1(t_c, x_c(:,7),   tp, 'linear');
    q_p = interp1(t_c, x_c(:,8:11),tp, 'linear')';
    % Pad biases and mass into 17-state vector
    x_hat_perfect = [r_p; v_p; q_p; zeros(6,1); m_p];

    [T_c, n_c, te] = guidance_track_reference(tp, x_hat_perfect, reference, params);
    fprintf('t=%5.1fs  T=%7.1f N  n=[%+.3f %+.3f %+.3f]  eR=%.2fm  eV=%.3fm/s\n', ...
        tp, T_c, n_c(1), n_c(2), n_c(3), te.pos_err, te.vel_err);
end
fprintf('  (Errors should be ~0 when fed the perfect reference state)\n');

%%  SECTION 4: Actuator saturation check 
fprintf('\n--- SECTION 4: Actuator Saturation ---\n');

if exist('plant_u_log','var')
    t_ulog = plant_u_log(:,1);
    T_log  = plant_u_log(:,2);
    a_log  = plant_u_log(:,3) * 180/pi;   % rad → deg
    b_log  = plant_u_log(:,4) * 180/pi;

    lim_a  = params.alpha_max * 180/pi;
    lim_b  = params.beta_max  * 180/pi;

    pct_a_sat = 100 * mean(abs(a_log) >= lim_a * 0.99);
    pct_b_sat = 100 * mean(abs(b_log) >= lim_b * 0.99);
    fprintf('Alpha saturation: %.1f%% of logged steps\n', pct_a_sat);
    fprintf('Beta  saturation: %.1f%% of logged steps\n', pct_b_sat);

    if pct_a_sat > 20 || pct_b_sat > 20
        fprintf('  *** WARNING: Gimbal is saturating heavily — attitude cannot follow guidance ***\n');
        fprintf('  FIX CANDIDATES:\n');
        fprintf('    1. Reduce LQI Q(1,1) (angle gain) in params.lqi.Q\n');
        fprintf('    2. Increase LQI R (control cost) in params.lqi.R\n');
        fprintf('    3. Increase alpha_max / beta_max if hardware allows\n');
    end

    figure('Name','Actuator Commands','NumberTitle','off');
    subplot(3,1,1); plot(t_ulog, T_log);
    ylabel('T (N)'); title('Thrust to Plant'); grid on;
    yline(params.Tmin,'r--'); yline(params.Tmax,'r--');

    subplot(3,1,2); plot(t_ulog, a_log);
    ylabel('\alpha (deg)'); title('Gimbal Pitch to Plant'); grid on;
    yline([ lim_a, -lim_a],'r--');

    subplot(3,1,3); plot(t_ulog, b_log);
    ylabel('\beta (deg)');  title('Gimbal Roll to Plant');  grid on;
    yline([ lim_b, -lim_b],'r--');
    xlabel('Time (s)');
else
    fprintf('plant_u_log not found — simulate first, or check plant_dynamics.m debug block.\n');
end

%% SECTION 5: LQI gain sanity check
fprintf('\n--- SECTION 5: LQI Gain Schedule ---\n');

K_t = params.K_table;
fprintf('%-8s  %-10s  %-8s  %-8s\n','Mass(kg)','K_theta','K_q','K_i');
for k = 1:size(K_t.K_pitch,1)
    fprintf('%-8.2f  %-10.4f  %-8.4f  %-8.4f\n', ...
        K_t.m_grid(k), K_t.K_pitch(k,1), K_t.K_pitch(k,2), K_t.K_pitch(k,3));
end

% Compute open-loop bandwidth at nominal mass
m_nom   = params.m0;
Ixy_nom = params.Ixy_per_kg * m_nom;
T_nom   = m_nom * params.g0;
b_nom   = T_nom * params.ls_nominal / Ixy_nom;
K_nom   = interp1(K_t.m_grid, K_t.K_pitch, m_nom, 'linear');
fprintf('\nAt m=%.2f kg:  b_pitch=%.3f rad/s² per rad\n', m_nom, b_nom);
fprintf('K = [%.3f, %.3f, %.3f]\n', K_nom(1), K_nom(2), K_nom(3));

% Check closed-loop poles
A_cl = [0, 1, 0; -K_nom(1)*b_nom, -K_nom(2)*b_nom, -K_nom(3)*b_nom; 1, 0, 0];
poles = eig(A_cl);
fprintf('Closed-loop poles: ');
fprintf('%.2f+%.2fj  ', [real(poles), imag(poles)]');
fprintf('\n');
if any(real(poles) >= 0)
    fprintf('  *** ERROR: UNSTABLE POLES — attitude loop will diverge ***\n');
else
    fprintf('  Poles OK (all LHP). Slowest time constant: %.2f s\n', ...
        -1/max(real(poles)));
end

%% SECTION 6: Mode A PD gain bandwidth check 
fprintf('\n--- SECTION 6: Mode A PD Gain Check ---\n');
Kp = params.modeA.Kp;
Kd = params.modeA.Kd;
fprintf('Kp = %.4f s^-2,  Kd = %.4f s^-1\n', Kp, Kd);

% Natural frequency and damping of the outer position loop
wn_pos  = sqrt(Kp);
zeta_pos = Kd / (2 * wn_pos);
fprintf('Outer loop: wn=%.3f rad/s (%.3f Hz),  zeta=%.3f\n', ...
    wn_pos, wn_pos/(2*pi), zeta_pos);

if wn_pos < 0.1
    fprintf('  *** WARNING: Kp very small — outer loop very slow (settling > %.0f s) ***\n', ...
        4/(zeta_pos*wn_pos));
    fprintf('  Try: params.modeA.Kp = 0.5; params.modeA.Kd = 1.2;\n');
end
if zeta_pos < 0.5
    fprintf('  *** WARNING: Under-damped outer loop — oscillations likely ***\n');
    fprintf('  Increase Kd to at least %.2f\n', 2*0.7*wn_pos);
end

%% SECTION 7: Closed-loop position tracking (from logsout) 
fprintf('\n SECTION 7: Tracking Error from Simulation Logs \n');
% This section expects workspace variables written by a Simulink To-Workspace
% block (or run_nominal.m). Adjust variable names to match your logging setup.

% Common variable names — edit to match yours:
var_t     = 'tout';          % time vector
var_r_hat = 'x_hat_log';     % EKF state log  (Nx17)
var_r_true= 'x_true_log';    % plant state log (Nx14)

if evalin('base', sprintf('exist(''%s'',''var'')', var_t))
    tout      = evalin('base', var_t);
    x_hat_log = evalin('base', var_r_hat);
    x_true_log= evalin('base', var_r_true);

    r_hat = x_hat_log(:,1:3);
    r_true= x_true_log(:,1:3);

    % Interpolate reference to simulation time grid
    r_ref_interp = interp1(t_c, x_c(:,1:3), tout, 'linear', 'extrap');

    pos_err_ekf  = vecnorm(r_ref_interp - r_hat,  2, 2);
    pos_err_true = vecnorm(r_ref_interp - r_true, 2, 2);

    figure('Name','Tracking Errors','NumberTitle','off');
    subplot(2,1,1);
    plot(tout, pos_err_ekf,  'b', tout, pos_err_true, 'r--');
    legend('EKF vs ref','Truth vs ref');
    ylabel('Position error (m)'); title('Tracking Error');
    yline(5,'k--','5m limit'); grid on;

    subplot(2,1,2);
    plot(tout, r_hat(:,3), 'b', tout, r_true(:,3), 'r--', ...
         interp1(t_c, t_c, tout, 'linear'), r_ref_interp(:,3), 'k:');
    legend('EKF Z','True Z','Ref Z');
    ylabel('Altitude (m)'); xlabel('Time (s)'); grid on;

    fprintf('Max position error (EKF vs ref): %.2f m\n',  max(pos_err_ekf));
    fprintf('Max position error (true vs ref): %.2f m\n', max(pos_err_true));

    if max(pos_err_true) > max(pos_err_ekf) * 2
        fprintf('  *** EKF tracking fine but PLANT diverges — check gains/actuators ***\n');
    elseif max(pos_err_ekf) > 5
        fprintf('  *** EKF error large — check GPS dropout rate and baro noise ***\n');
    end
else
    fprintf('Simulation log variables not found (%s, %s, %s).\n', ...
        var_t, var_r_hat, var_r_true);
    fprintf('Add To-Workspace blocks in Simulink and re-run, or adapt variable names above.\n');
end

%%  SECTION 8: n_hat_cmd direction alignment check 
fprintf('\n--- SECTION 8: Thrust Direction vs Attitude Consistency ---\n');
% Does attitude_reference invert n_hat_cmd consistently with the plant?
t_test   = 10;   % pick a time mid-ascent
r_pt     = interp1(t_c, x_c(:,1:3), t_test, 'linear')';
v_pt     = interp1(t_c, x_c(:,4:6), t_test, 'linear')';
m_pt     = interp1(t_c, x_c(:,7),   t_test, 'linear');
q_pt     = interp1(t_c, x_c(:,8:11),t_test, 'linear')';
x_hat_t  = [r_pt; v_pt; q_pt; zeros(6,1); m_pt];

[~, n_hat] = guidance_track_reference(t_test, x_hat_t, reference, params);
[theta_r, phi_r] = attitude_reference(n_hat);

fprintf('At t=%.1f s:\n', t_test);
fprintf('  n_hat_cmd = [%.4f, %.4f, %.4f]  (should be unit: norm=%.6f)\n', ...
    n_hat(1), n_hat(2), n_hat(3), norm(n_hat));
fprintf('  theta_ref = %.3f deg,  phi_ref = %.3f deg\n', ...
    theta_r*180/pi, phi_r*180/pi);

% Reconstruct the body-z axis from the demanded attitude and check it
% agrees with n_hat (they describe the same geometric direction)
q_test = [cos(theta_r/2)*cos(phi_r/2);
          sin(phi_r/2)*cos(theta_r/2);
          cos(phi_r/2)*sin(theta_r/2);
          sin(theta_r/2)*sin(phi_r/2)];
q_test = q_test / norm(q_test);
R_bi_t = quat2rotm(q_test');           % MATLAB toolbox
z_body_in_inertial = R_bi_t * [0;0;1]; % body-z in inertial frame

dot_product = dot(n_hat, z_body_in_inertial);
fprintf('  Alignment (dot product, should be ~1.0): %.4f\n', dot_product);
if dot_product < 0.95
    fprintf('  *** WARNING: attitude_reference decomposition is inconsistent with n_hat ***\n');
    fprintf('  The sign convention in attitude_reference.m may not match plant_dynamics.m\n');
end

fprintf('\n--- DIAGNOSIS COMPLETE ---\n');
fprintf('Next steps based on results above:\n');
fprintf('  Saturating gimbal  → retune LQI (Section 5 advice)\n');
fprintf('  Outer loop too slow → increase Kp/Kd (Section 6 advice)\n');
fprintf('  Reference bad      → re-run ICLOCS and check bounds\n');
fprintf('  Alignment mismatch → check sign convention in attitude_reference.m\n');



%  CONFIRM_BUGS.M  — paste into command window
% ================================================================

%% Bug 1: Where does T_ref drop below Tmin?
[t_c, ia] = unique(reference.t, 'stable');
u_c  = reference.u(ia,:);
x_c  = reference.x(ia,:);
ph_c = reference.phase(ia);

T_ref_all = u_c(:,1);
below_Tmin = T_ref_all < params.Tmin;

fprintf('=== BUG 1: T_ref < Tmin ===\n');
fprintf('Fraction of reference with T_ref < Tmin: %.1f%%\n', ...
    100*mean(below_Tmin));
fprintf('T_ref at t=0 : %.2f N  (Tmin = %.0f N)\n', ...
    interp1(t_c, T_ref_all, 0, 'linear'), params.Tmin);
fprintf('Mean T_ref during coast (phases 1-3): %.2f N\n', ...
    mean(T_ref_all(ph_c < 4)));

% Plot to see coast windows clearly
figure('Name','Bug 1 — Coast Windows','NumberTitle','off');
subplot(2,1,1);
plot(t_c, T_ref_all, 'b'); hold on;
yline(params.Tmin, 'r--', 'Tmin=32 N', 'LabelHorizontalAlignment','left');
fill([t_c(below_Tmin); flipud(t_c(below_Tmin))], ...
    [zeros(sum(below_Tmin),1); params.Tmin*ones(sum(below_Tmin),1)], ...
    'r', 'FaceAlpha', 0.15, 'EdgeColor','none');
ylabel('T_{ref} (N)'); title('Reference Thrust — Red shading = T_{ref} < Tmin (WRONG)');
grid on; legend('T_{ref}','Tmin','Coast zone');

subplot(2,1,2);
plot(t_c, double(ph_c));
ylabel('Phase'); xlabel('t (s)'); grid on;
yticks(1:4); ylim([0.5 4.5]);
title('Phase schedule');

%% Bug 2: What does guidance actually command vs what plant receives?
% Call guidance at several coast-phase times
fprintf('\n=== BUG 2: Guidance command vs actuator output ===\n');
fprintf('%-6s  %-6s  %-10s  %-10s  %-10s\n', ...
    'Time','Phase','T_ref(N)','T_guidance','After clamp');

t_check = [0, 1, 5, 10, 20, 30, 35];
for tp = t_check
    if tp > t_c(end), continue; end
    r_p = interp1(t_c, x_c(:,1:3), tp)';
    v_p = interp1(t_c, x_c(:,4:6), tp)';
    m_p = interp1(t_c, x_c(:,7),   tp);
    q_p = interp1(t_c, x_c(:,8:11),tp)';
    ph  = interp1(t_c, double(ph_c), tp, 'nearest');
    Tr  = interp1(t_c, T_ref_all, tp);
    x_hat_p = [r_p; v_p; q_p; zeros(6,1); m_p];
    [T_g, ~] = guidance_track_reference(tp, x_hat_p, reference, params);
    T_after  = max(params.Tmin, min(params.Tmax, T_g));
    fprintf('t=%4.0fs  ph=%d  T_ref=%6.1f  T_guid=%7.1f  T_plant=%7.1f  %s\n', ...
        tp, ph, Tr, T_g, T_after, ...
        iff(T_after > Tr*2, '<-- FORCED UP by clamp', ''));
end

%% Bug 3: What attitude does guidance demand at t=0?
fprintf('\n=== BUG 3: Initial attitude demand ===\n');
x0_perfect = [0;0;0; 0;0;0; [1;0;0;0]; zeros(6,1); params.m0];
[~, n0] = guidance_track_reference(0, x0_perfect, reference, params);
[theta0, phi0] = attitude_reference(n0);
fprintf('At t=0, guidance demands:\n');
fprintf('  n_hat_cmd = [%.3f, %.3f, %.3f]\n', n0(1), n0(2), n0(3));
fprintf('  theta_ref = %.2f deg  (attitude loop error at t=0: %.2f deg)\n', ...
    theta0*180/pi, theta0*180/pi);  % vehicle starts at theta=0
fprintf('  phi_ref   = %.2f deg  (attitude loop error at t=0: %.2f deg)\n', ...
    phi0*180/pi, phi0*180/pi);
fprintf('  alpha_cmd (before sat) = %.2f deg  (limit = ±%.1f deg)\n', ...
    -params.K_table.K_pitch(1,1)*theta0*180/pi, params.alpha_max*180/pi);

if abs(-params.K_table.K_pitch(1,1)*theta0) > params.alpha_max
    fprintf('  *** Gimbal IMMEDIATELY saturates at t=0 ***\n');
end