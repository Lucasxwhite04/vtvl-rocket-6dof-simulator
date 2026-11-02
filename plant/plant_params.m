function params = plant_params()
% PLANT_PARAMS  Returns the nominal vehicle parameters struct for the VTVL prototype.
%
% Called once at startup; the returned struct is passed explicitly to every
% function that needs vehicle data 
%
% ---- Vehicle summary ---
%   Body:           1.2 m × 0.15 m cylinder
%   Dry mass:       8 kg
%   Propellant:     8.88 kg  →  initial total mass 16.88 kg
%   Max thrust:     3200 N,  min thrust 32 N (1% deep-throttle)
%   TVC gimbal:     2-axis (pitch α, roll β), ±5°, no yaw authority
%
% ---- State vector convention (14×1 column) ---
%   x = [X; Y; Z; vx; vy; vz; m; q0; q1; q2; q3; p; q_rate; r_rate]
%
%   Index  Symbol    Description                     Units
%   -----  ------    -----------                     -----
%    1      X        Inertial position X             m
%    2      Y        Inertial position Y             m
%    3      Z        Inertial position Z  (Z-up)     m
%    4      vx       Inertial velocity X             m/s
%    5      vy       Inertial velocity Y             m/s
%    6      vz       Inertial velocity Z             m/s
%    7      m        Total mass                      kg
%    8      q0       Quaternion scalar               —
%    9      q1       Quaternion vector-x             —
%   10      q2       Quaternion vector-y             —
%   11      q3       Quaternion vector-z             —
%   12      p        Body roll  rate (about body-x)  rad/s
%   13      q_rate   Body pitch rate (about body-y)  rad/s
%   14      r_rate   Body yaw   rate (about body-z)  rad/s  [free — yaw uncontrolled]
%
%   Quaternion convention: q = [q0; q1; q2; q3], q0 scalar-first (Hamilton).
%   Angles are in RADIANS internally; convert to degrees only at plot/UI boundaries.
%
% ---- Control vector (3×1) ---
%   u = [T; alpha; beta]
%   T     - thrust magnitude (N),        bounds [Tmin, Tmax]
%   alpha - gimbal pitch angle (rad),    bounds ±alpha_max
%   beta  - gimbal roll  angle (rad),    bounds ±beta_max
%
% ---- Outputs ---
%   params  struct (all SI units unless explicitly noted)
% L. White, University of Sheffield, May 2026




%%  GEOMETRY
params.length   = 1.2;                        % vehicle length   (m)
params.diameter = 0.15;                       % vehicle diameter (m)
params.radius   = params.diameter / 2;        % vehicle radius   (m)



%%  MASS PROPERTIES

params.ms = 8.00;                              % dry mass (structural + avionics) (kg)
params.mp = 8.88;                              % propellant mass (kg)
params.m0 = params.ms + params.mp;             % initial total mass (kg)  = 16.88

% --- Moments of inertia (thin uniform cylinder) ---
% Transverse (body-x = body-y, symmetric rocket):
%   I_xy = (1/12)*m*L^2  +  (1/4)*m*R^2
% Axial (body-z, spin axis):
%   I_z  = (1/2)*m*R^2
%
% Stored as per-unit-mass so plant_dynamics.m can recompute at each
% integration step as mass depletes:
%   Ixy(m) = params.Ixy_per_kg * m
%   Iz(m)  = params.Iz_per_kg  * m
%
% % TBD: update once avionics bay layout and tank geometry are finalised.
params.Ixy_per_kg = (1/12)*params.length^2 + (1/4)*params.radius^2;  % m²/kg
params.Iz_per_kg  = (1/2) *params.radius^2;                           % m²/kg

% --- IMU lever-arm ---
% Position of the IMU in the body frame (m).
% Used for lever-arm correction in imu_model.m and ekf_predict.m.
% Currently set to the body-frame origin (CoM) — no lever-arm effect.
% % TBD: update to avionics bay offset [0; 0; z_imu] once bay is located.
params.r_imu = [0; 0; 0];    % (3×1) body frame (m)


%%  PROPULSION
params.Tmin      = 32;                          % minimum thrust (N)  — 1% of Tmax
params.Tmax      = 3200;                        % maximum thrust (N)
params.alpha_max =  5 * pi/180;                 % gimbal pitch limit (rad)
params.beta_max  =  5 * pi/180;                 % gimbal roll  limit (rad)

% Effective exhaust velocity  c = Isp * g0
% Bipropellant Isp ≈ 270 s (conservative prototype estimate).
params.Isp = 270;                               % specific impulse (s)
params.g0  = 9.80665;                           % standard gravity (m/s²)
params.c   = 1488;                              % effective exhaust velocity (m/s)

% TVC moment arm nominal value.
% compute_moment_arm() in plant_dynamics.m interpolates between ls_nominal
% (fully loaded) and ls_nominal + 0.05 m (dry) as propellant burns.
% CoM ≈ mid-body for a uniform vehicle → ls ≈ length/2.
% % TBD: update with propellant tank geometry when finalised.
params.ls_nominal = params.length / 2;          % (m)


%%  AERODYNAMICS
% Drag model:  D_i = -0.5 * rho * CD_i * A_i * V_i * |V_i|
% Smooth C¹ form — no singularity at V = 0 (Build Plan §5.2).
%
% Lateral axes (body-x, body-y):  CDxy, A_lateral
% Axial axis   (body-z):          CDz,  A_axial
params.CDxy      = 0.8;                         % lateral drag coefficient  (—)
params.CDz       = 0.3;                         % axial   drag coefficient  (—)
params.A_lateral = params.length * params.diameter;   % projected lateral area (m²)
params.A_axial   = pi * params.radius^2;              % circular cross-section (m²)

% ISA sea-level density (treated as constant: mission altitude < 300 m,
% Δρ/ρ < 3.5%, negligible for this vehicle).
params.rho = 1.225;                             % air density (kg/m³)


%%  GRAVITY
% Constant-gravity approximation valid for h < 300 m.
% Z-up inertial convention: g acts in the −Z direction.
params.grav = [0; 0; -params.g0];              % (3×1) m/s²

%%  SENSOR AND EKF PARAMETERS

% --- GPS ---
params.gps.pos_noise_std      = [0.5; 0.5; 0.5];    % position noise 1-sigma [m]
params.gps.vel_noise_std      = [0.05; 0.05; 0.05];  % velocity noise 1-sigma [m/s]
params.gps.dropout_prob       = 0.10;                % 10% random dropout rate
params.gps.rate_hz            = 10;                  % GPS update rate        [Hz]
params.gps.mission_duration_s = 65;                  % pre-allocation window  [s]

% --- IMU ---
params.imu.accel_bias        = [0; 0; 0];        % mean turn-on bias  [m/s²]
params.imu.accel_bias_sigma  = [0; 0; 0];        % zeroed for diagnostic
params.imu.accel_noise_std   = [0.01; 0.01; 0.01];
params.imu.accel_scale_err   = [0; 0; 0];
params.imu.gyro_bias         = [0; 0; 0];        % mean turn-on bias  [rad/s]
params.imu.gyro_bias_sigma   = [0; 0; 0];        % zeroed for diagnostic
params.imu.gyro_noise_std    = [0.001; 0.001; 0.001];
params.imu.gyro_scale_err    = [0; 0; 0];
% --- Barometer ---
params.baro.noise_std  = 0.2;    % per-sample white noise 1-sigma  [m]
params.baro.bias       = 0;
params.baro.bias_sigma = 0;   % flight-to-flight bias variation  [m]
                                 % bias_rt = bias + bias_sigma*randn()
                                 % gives realistic run-to-run scatter
% --- EKF tuning ---
params.ekf.mass_init_std          = 0.888;     % initial mass uncertainty [kg] (~10%)
params.ekf.mass_process_noise_std = 0.01;      % mass random walk noise   [kg/s^0.5]
params.ekf.accel_noise_std        = [0.01; 0.01; 0.01]; % [m/s²]
params.ekf.gyro_noise_std         = [0.001; 0.001; 0.001]; % [rad/s]
params.ekf.gyro_bias_rw           = [1e-4; 1e-4; 1e-4];   % bias random walk [rad/s/s^0.5]
params.ekf.accel_bias_rw          = [1e-3; 1e-3; 1e-3];   % bias random walk [m/s²/s^0.5]

%%  MISSION GEOMETRY
params.r_pad    = [200; -200; 0];               % landing pad position (m), inertial
params.n_phases = 4;                            % active mission phases (3 default, 4 with hover)
                                                % set to 4 and re-enable HOVER_DISABLED blocks
                                                % to switch to the 4-phase mission


%%  CLOSED-LOOP GNC PARAMETERS

% --- Guidance / MPC ---
% Pipeline latency: time between MEKF output and first actuator command.
% The State Predictor (guidance/state_predictor.m) projects x_hat forward
% by this delay before passing it to mpc_solve.m (G-FOLD pattern, §12.14).
% Measure on the target hardware and update; default 100 ms for desktop.
params.t_plan_delay = 0.100;                    % planning pipeline delay (s)

% --- MODE A GUIDANCE GAINS ---  
% PD correction gains for reference-tracking guidance.
% Start conservative; increase gradually once nominal loop is verified.
% Kp [s^-2]: position error gain — increase if vehicle drifts slowly.
% Kd [s^-1]: velocity error gain — dominates short-term damping.
params.modeA.Kp = 0.3;      % position error gain  [s^-1]
params.modeA.Kd = 0.6;      % velocity error gain  [s^-1]
% --- Yaw monitor ---
% yaw_monitor.m (control/) raises an alarm flag when |ψ| exceeds this
% threshold. Does NOT attempt yaw control — the 2-axis gimbal has no
% body-z torque authority.
params.yaw_alarm_threshold = 10 * pi/180;       % 10° alarm threshold (rad)

% --- Non-dimensionalisation constant ---
% Used inside mpc_solve.m and scp_iteration.m to condition the SOCP.
% Re-dimensionalise solver outputs before returning to the inner loop.
params.U_L = 100;                               % length scale    (m)
params.U_T = sqrt(params.U_L / params.g0);      % time scale      (s)  ≈ 3.19 s
params.U_M = params.m0;                         % mass scale      (kg) = 16.88

%%  LQI ATTITUDE LOOP TUNING WEIGHTS 
% Used by lqi_gains.m. Tune these, not the gains directly.
% Q(1): angle error weight  — bigger = faster response, more gimbal use
% Q(2): rate weight         — bigger = more damping, less overshoot  
% Q(3): integral weight     — bigger = faster steady-state correction
% R:    control effort      — bigger = gentler gimballing
params.lqi.Q = diag([100, 100, 1]);
params.lqi.R = 1000;


%%  MODE B GUIDANCE — MPC parameters 


% --- Guidance loop rate ---
% Mode B replans every dt_guidance seconds. Start at 1 Hz; tighten
% to 2 Hz (0.5 s) once nominal Step 7 passes.
params.dt_guidance = 1.0;                       % MPC replan period (s)

% --- SCP outer loop ---
params.scp.N              = 30;                 % discretisation nodes
params.scp.i_max          = 8;                  % max SCP iterations per call
params.scp.nu_tol         = 1e-4;               % virtual-control inf-norm tol [non-dim]
params.scp.delta_tol      = 1e-3;               % trust-region inf-norm tol    [non-dim]
params.scp.lcvx_tol       = 0.05;               % LCvx slack tol  |Γ - ||u|||  [N]
params.scp.max_rejections = 3;                  % consecutive rejects → Mode A fallback

% --- SOCP cost weights ----
params.scp.w_nu     = 1e3;                      % virtual-control penalty
params.scp.w_delta  = 1e2;                      % trust-region slack penalty
params.scp.w_t_desc = 0.02;                        % time penalty — descent (free t_f)
params.scp.w_t_land = 0;                        % time penalty — landing (fixed t_f)
params.scp.cost_tol = 1e-4;                     % fuel-cost change tolerance [kg] (~1 mg)



% --- Path constraints ---
% Glideslope half-angle (cone vertex at landing pad, opening upward).
% Start loose at 60°; tighten to 50° once Step 7 lands cleanly.
params.gamma_gs_deg = 60;                       % glideslope half-angle (deg)
params.gamma_gs     = params.gamma_gs_deg * pi/180;  % (rad)

% --- Terminal box tolerances (used by validity_check.m) ---
params.terminal_pos_tol = 2.5;                  % ‖r_f - r_target‖ (m)
params.terminal_vel_tol = 0.2;                  % ‖v_f - v_target‖ (m/s)

% --- Time bounds per phase ---
% Phase 3 (descent): free final time σ ∈ [5, 30] s
% Phase 4 (landing): fixed σ = 3.0 s  (from reference_trajectory.mat)
params.scp.sigma_desc_low  = 5.0;               % descent  lower bound (s)
params.scp.sigma_desc_high = 30.0;              % descent  upper bound (s)
params.scp.sigma_land      = 3.0;               % landing  fixed value (s)
params.scp.delta_max       = 1.0;               % trust-region radius (non-dim L2 norm)
params.scp.z_cap_margin = 1.0;                  % metres above start altitude allowed

% --- Derived non-dim quantities (for convenience) ----
% U_L, U_T, U_M already defined above. Add the dependent units:
params.U_F = params.U_M * params.U_L / params.U_T^2;  % force scale (N) ≈ 165.5
params.U_V = params.U_L / params.U_T;                  % velocity (m/s) ≈ 31.3
params.U_A = params.U_L / params.U_T^2;                % accel (m/s²) = g0

%%  WIND MODEL PARAMETERS  

% Nominal run: all zeros — no wind, no turbulence.
params.wind.enable              = false;  % master on/off switch
params.wind.mean_speed          = 0;      % mean wind speed          [m/s]
params.wind.direction_deg       = 270;    % wind FROM direction       [deg] West
params.wind.ref_altitude        = 10;     % reference altitude        [m]
params.wind.shear_exp           = 0.14;   % Hellmann power law exp    [-]
params.wind.turbulence_intensity= 0;      % gust intensity multiplier [-]
params.wind.gust_sigma          = 0;      % turbulence std dev        [m/s]
params.wind.gust_bandwidth      = 0.5;    % turbulence bandwidth      [rad/s]
params.wind.seed                = 42;     % RNG seed for repeatability


%%  SIMULATION SEED
params.sim_seed = 42;    % RNG seed for all sensor noise

%%  MONTE CARLO DISPERSIONS  

% ---- Uncertainty bounds (fractional, ±) ---------------------------
params.dispersions.bounds.mp_frac   = 0.10;    % ±10% propellant mass
params.dispersions.bounds.CDxy_frac = 0.20;    % ±20% lateral drag coefficient
params.dispersions.bounds.CDz_frac  = 0.20;    % ±20% axial   drag coefficient
params.dispersions.bounds.Isp_frac  = 0.02;    % ±2%  specific impulse / exhaust velocity

% ---- Per-run scale factors (nominal = 1.0, overwritten by run_montecarlo.m) --
% plant_dynamics.m applies these as:
%   rho_applied  = params.rho  * dispersions.scale.rho_scale
%   CDxy_applied = params.CDxy * dispersions.scale.CDxy_scale   etc.
params.dispersions.scale.rho_scale  = 1.0;     % air density      (—)
params.dispersions.scale.CDxy_scale = 1.0;     % lateral drag CD  (—)
params.dispersions.scale.CDz_scale  = 1.0;     % axial    drag CD (—)
params.dispersions.scale.c_scale    = 1.0;     % exhaust velocity (—)


% Load reference trajectory into params for Simulink access
ref_path = 'C:\Users\lucas\OneDrive\Documents\FYP\A4 - Rocket\V7\offline\reference_trajectory.mat';
if exist(ref_path, 'file')
    tmp = load(ref_path, 'reference');
    params.reference = tmp.reference;
end


%%  PRECOMPUTED GAINS  (computed once at startup)
params.K_table = lqi_gains(params);   % LQI gain schedule table

%% compile the SOCP optimiser once at startup (~1 s)
params.socp_optimizer = build_socp(params);

end