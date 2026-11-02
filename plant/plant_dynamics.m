function xdot = plant_dynamics(x, u, t, params, wind)
% PLANT_DYNAMICS  6-DoF equations of motion for the VTVL prototype.
%
% Quaternion-based truth model. Called by the Simulink plant block at 1000 Hz
% (ode45 variable-step). This is the unified ground-truth plant; it is NOT
% used directly by the MPC (which uses the 3-DoF translational model).
%
% --- State vector (14×1) ---
%
%   Index  Symbol     Description                         Units
%   -----  ------     -----------                         -----
%    1      X         Inertial position X                 m
%    2      Y         Inertial position Y                 m
%    3      Z         Inertial position Z  (Z-up)         m
%    4      vx        Inertial velocity X                 m/s
%    5      vy        Inertial velocity Y                 m/s
%    6      vz        Inertial velocity Z                 m/s
%    7      m         Total mass                          kg
%    8      q0        Quaternion scalar part              —
%    9      q1        Quaternion vector part x            —
%   10      q2        Quaternion vector part y            —
%   11      q3        Quaternion vector part z            —
%   12      p         Body roll  rate (about body-x)      rad/s
%   13      q_rate    Body pitch rate (about body-y)      rad/s
%   14      r_rate    Body yaw   rate (about body-z)      rad/s  
%
%
% --- Control vector (3×1) ---
%
%   Index  Symbol     Description                         Units
%   -----  ------     -----------                         -----
%    1      T         Thrust magnitude                    N      [32, 3200]
%    2      alpha     Gimbal pitch deflection             rad    [−5°, +5°]
%    3      beta      Gimbal roll  deflection             rad    [−5°, +5°]
%
% --- Inputs ---
%   x      (14×1)   State vector (see above)
%   u      (3×1)    Control inputs [T; alpha; beta]
%   t      (scalar) Current simulation time (s)  — retained for Simulink API
%   params           Struct from plant_params.m
%   wind   (3×1)    Wind velocity in inertial frame (m/s).
%                   Pass zeros(3,1) if no wind model is connected.
%
% --- Outputs ---
%   xdot   (14×1)   State time-derivative
%
% L. White, University of Sheffield, May 2026


%%  STEP 1 - Apply parameter dispersions (Monte Carlo support)
% Local copies of aerodynamic / propulsion parameters.
% If params.dispersions exists and is non-empty, scale factors are
% applied multiplicatively. This leaves params unmodified (no side-effects),
% so parfor workers do not interfere with each other.

rho  = params.rho;    % air density (kg/m³)
CDxy = params.CDxy;   % lateral drag coefficient (—)
CDz  = params.CDz;    % axial   drag coefficient (—)
c    = params.c;      % effective exhaust velocity (m/s)

if isfield(params, 'dispersions') && ~isempty(params.dispersions)
    d = params.dispersions;
    % Each field is optional — only applied if present.
    if isfield(d, 'rho_scale'),  rho  = rho  * d.rho_scale;  end  % e.g. ±20%
    if isfield(d, 'CDxy_scale'), CDxy = CDxy * d.CDxy_scale; end  % e.g. ±20%
    if isfield(d, 'CDz_scale'),  CDz  = CDz  * d.CDz_scale;  end  % e.g. ±20%
    if isfield(d, 'c_scale'),    c    = c    * d.c_scale;    end  % propulsion Isp
end


%%  STEP 1 - Unpack and validate state


% Position (inertial frame, Z-up)
X  = x(1);
Y  = x(2);
Z  = x(3);

% Velocity (inertial frame)
vx = x(4);
vy = x(5);
vz = x(6);

% Total mass
m  = x(7);

% --- Quaternion attitude — renormalise immediately ---
% Integrator drift causes the norm to deviate from 1 over thousands
% of steps. Renormalising here costs negligible computation and keeps
% the rotation matrix orthonormal throughout the simulation.
q0 = x(8);
q1 = x(9);
q2 = x(10);
q3 = x(11);

q_norm = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
if q_norm > eps   % guard against degenerate zero-quaternion
    q0 = q0 / q_norm;
    q1 = q1 / q_norm;
    q2 = q2 / q_norm;
    q3 = q3 / q_norm;
end

% --- Body-frame angular rates ---
p      = x(12);   % roll  rate about body-x (rad/s)
q_rate = x(13);   % pitch rate about body-y (rad/s)
r_rate = x(14);   % yaw   rate about body-z (rad/s)  — free to evolve, not controlled


%%  STEP 2 — Clamp controls to physical limits
% Safety clamp: protect dynamics from optimiser outputs that violate
% actuator bounds. Applied before any physics computation.

T     = min(max(u(1), 0),               params.Tmax);        % thrust (N)
alpha = min(max(u(2), -params.alpha_max), params.alpha_max);  % gimbal pitch (rad)
beta  = min(max(u(3), -params.beta_max),  params.beta_max);   % gimbal roll  (rad)

% (TEMP DEBUG LINES)
persistent plant_u_log_
if isempty(plant_u_log_), plant_u_log_ = []; end
if size(plant_u_log_, 1) < 500   % log first 500 calls only
    plant_u_log_(end+1, :) = [t, u(:)'];   % columns: [time, T, alpha, beta]
    assignin('base', 'plant_u_log', plant_u_log_);
end
% Once the propellant is exhausted, thrust cannot be commanded.
% This guards against numerical overshoot near m = ms.
if m <= params.ms
    T = 0;
end


%%  STEP 3 — Rotation matrix: body → inertial  (R_bi)
% R_bi maps a vector expressed in the rocket body frame to the
% inertial (Earth) frame.  Derived from unit quaternion q = [q0;q1;q2;q3]
% (scalar-first, Hamilton convention).

R_bi = [ 1 - 2*(q2^2 + q3^2),    2*(q1*q2 - q0*q3),    2*(q1*q3 + q0*q2) ;
              2*(q1*q2 + q0*q3),  1 - 2*(q1^2 + q3^2),  2*(q2*q3 - q0*q1) ;
              2*(q1*q3 - q0*q2),  2*(q2*q3 + q0*q1),    1 - 2*(q1^2 + q2^2) ];

% Inverse rotation (inertial → body):
R_ib = R_bi';   % (R_bi is orthonormal, so transpose = inverse)



%%  STEP 4 — Thrust vector in inertial frame
% The gimballed nozzle deflects the thrust direction relative to the
% rocket body axis.  Sign convention:
%   +alpha tilts toward +body-x → generates +pitch moment about body-y (My > 0)
%   +beta  tilts toward +body-y → generates +roll  moment about body-x (Mx < 0)
% The dominant component is axial (along +body-z).
%
%   Tr_body = [ T·sin(alpha) ; T·sin(beta) ; T·cos(alpha)·cos(beta) ]
%
% Body-x convention: when the rocket is upright, body-x points in the
% −X_inertial direction, so the x-component carries a sign flip before
% rotation to the inertial frame.

Tr_body = [ T * sin(alpha)          ;   % lateral-x (pitch) component
            T * sin(beta)           ;   % lateral-y (roll)  component
            T * cos(alpha)*cos(beta) ]; % axial      (dominant)      component

% Apply sign flip on body-x component over from ICLOCS offline model for consistency at the interfaces.
Tr_body_signed = [-Tr_body(1); Tr_body(2); Tr_body(3)];

% Rotate from body frame to inertial frame
Te = R_bi * Tr_body_signed;   % (3×1) thrust force in inertial frame (N)


%%  STEP 5 — Aerodynamic drag (smooth C¹ formulation)
% Smooth drag: D_i = −0.5 · ρ · CD_i · A_i · V_i · |V_i|
%
% This replaces the non-smooth v./abs(v) formulation from the offline
% ICLOCS phase dynamics.  The two forms are
% numerically identical for V_i ≠ 0; the smooth form avoids a gradient
% kink at v = 0 that causes IPOPT to slow down and ode45 to overshoot.
%
% The drag calculation is performed in the body frame because the
% aerodynamic coefficients CDxy and CDz differ by axis (asymmetric
% cross-section).  Result is rotated back to inertial frame.

% Relative airspeed: subtract wind to get the aerodynamically significant
% velocity (inertial frame first, then transform to body).
v_rel_inertial = [vx; vy; vz] - wind;   % (3×1) m/s, inertial frame

% Transform to body frame for component-wise drag computation
V_body = R_ib * v_rel_inertial;   % (3×1) m/s, body frame

% Per-axis drag in body frame.
% CDxy applies to lateral axes (body-x, body-y); CDz applies to the axial
% direction (body-z). A_lateral and A_axial are the projected cross-sections.
Dr_body = [ -0.5 * rho * CDxy * params.A_lateral * V_body(1) * abs(V_body(1)) ;
            -0.5 * rho * CDxy * params.A_lateral * V_body(2) * abs(V_body(2)) ;
            -0.5 * rho * CDz  * params.A_axial   * V_body(3) * abs(V_body(3)) ];

% Rotate drag back to inertial frame
De = R_bi * Dr_body;   % (3×1) drag force in inertial frame (N)


%%  STEP 6 — Translational dynamics  (Newton's 2nd law)
% ṙ = v
% v̇ = (Te + De) / m  +  g
%
% params.grav = [0; 0; -g0]  (Z-up inertial, g0 ≈ 9.806 m/s²)

a_inertial = (Te + De) / m + params.grav;   % (3×1) m/s²


%%  STEP 7 — Mass flow rate
% Tsiolkovsky propellant consumption:
%   ṁ = −T / c
%
% c = effective exhaust velocity (m/s). 

mdot = -T / c;   % kg/s  (negative: mass decreases with thrust)


%%  STEP 8 — TVC moments in body frame
%
% The gimballed nozzle displaces the thrust line by moment arm ls
% (CoM to nozzle exit).  ls changes as propellant burns and CoM shifts.

ls = compute_moment_arm(m, params);   % m   (see sub-function below)

% Roll  moment (about body-x):  Mx = −T · sin(β) · ls
%   (positive β tilts thrust toward +body-y, creating a clockwise roll
%    viewed from nose, negative Mx by right-hand rule)
Mx_tvc = -T * sin(beta)  * ls;   % N·m

% Pitch moment (about body-y):  My = +T · sin(α) · ls
%   (positive α tilts thrust toward +body-x, nose pitches up, My > 0)
My_tvc =  T * sin(alpha) * ls;   % N·m

% Yaw moment from TVC:  Mz_tvc = 0
%   The centerline 2-axis gimbal produces NO torque about body-z.
Mz_tvc = 0;

% --- Disturbance yaw moment ---
% Placeholder for future wind-drag asymmetry model.
% At present Mz_dist = 0; add an expression here once wind_model.m
% produces spatially asymmetric forces (e.g. gust profile across the
% vehicle length).
Mz_dist = 0;   % N·m

% Total moments
Mx = Mx_tvc;             % N·m
My = My_tvc;             % N·m
Mz = Mz_tvc + Mz_dist;  % N·m  (= 0 until disturbance model is added)


%%  STEP 9 — Angular dynamics  (Euler's equations, full 3-DoF)
% Euler's equations for a rigid body under applied moments:
%
%   Ix · ṗ       = Mx + (Iy − Iz) · q_rate · r_rate
%   Iy · q̇_rate = My + (Iz − Ix) · r_rate · p
%   Iz · ṙ_rate  = Mz + (Ix − Iy) · p · q_rate
%
% For this symmetric vehicle  Ix = Iy = Ixy  (transverse inertia):
%   ṗ       = [Mx + (Ixy − Iz) · q_rate · r_rate] / Ixy
%   q̇_rate = [My + (Iz − Ixy) · r_rate · p      ] / Ixy
%   ṙ_rate  = [Mz + 0                             ] / Iz
%            = Mz / Iz
Ixy = params.Ixy_per_kg * m;   % kg·m²  (pitch = roll, symmetric rocket)
Iz  = params.Iz_per_kg  * m;   % kg·m²  (axial spin — small for slender body)

pdot      = (Mx + (Ixy - Iz) * q_rate * r_rate) / Ixy;  % roll  angular accel (rad/s²)
qdot_rate = (My + (Iz - Ixy) * r_rate * p) / Ixy;  % pitch angular accel (rad/s²)
rdot_rate = (Mz) / Iz;   % yaw   angular accel (rad/s²)
%                                                          % cross-coupling = 0  (Ix = Iy)


%%  STEP 10 — Quaternion kinematics
% q̇ = ½ · Ξ(q) · ω_body
%
% Ξ(q) is the 4×3 kinematic matrix (scalar-first quaternion convention,
% Hamilton product form).  Derivation: Wie 2008, Eq. 4.4-12.
%
%        [ −q1   −q2   −q3 ]
% Ξ(q) = [  q0   −q3    q2 ]
%        [  q3    q0   −q1 ]
%        [ −q2    q1    q0 ]
%
% Verification (spot check):
%   q̇0 = ½(−q1·p − q2·q_rate − q3·r_rate)
%   q̇1 = ½( q0·p − q3·q_rate + q2·r_rate)
%   q̇2 = ½( q3·p + q0·q_rate − q1·r_rate)
%   q̇3 = ½(−q2·p + q1·q_rate + q0·r_rate)

omega_body = [p; q_rate; r_rate];   % (3×1) body angular velocity (rad/s)

Xi_q = [ -q1,  -q2,  -q3 ;
          q0,  -q3,   q2 ;
          q3,   q0,  -q1 ;
         -q2,   q1,   q0 ];          % (4×3) kinematic matrix

qdot_vec = 0.5 * Xi_q * omega_body; % (4×1) quaternion time-derivative


%%  STEP 11 — Assemble state derivative
xdot = zeros(14, 1);

% Translational kinematics
xdot(1)  = vx;               % Ẋ          (m/s)
xdot(2)  = vy;               % Ẏ          (m/s)
xdot(3)  = vz;               % Ż          (m/s)

% Translational dynamics
xdot(4)  = a_inertial(1);    % v̇x         (m/s²)
xdot(5)  = a_inertial(2);    % v̇y         (m/s²)
xdot(6)  = a_inertial(3);    % v̇z         (m/s²)

% Mass depletion
xdot(7)  = mdot;             % ṁ           (kg/s)

% Quaternion kinematics
xdot(8)  = qdot_vec(1);      % q̇0          (1/s)
xdot(9)  = qdot_vec(2);      % q̇1          (1/s)
xdot(10) = qdot_vec(3);      % q̇2          (1/s)
xdot(11) = qdot_vec(4);      % q̇3          (1/s)

% Angular dynamics
xdot(12) = pdot;             % ṗ  (roll  angular accel) (rad/s²)
xdot(13) = qdot_rate;        % q̈  (pitch angular accel) (rad/s²)
xdot(14) = rdot_rate;        % ṙ  (yaw   angular accel) (rad/s²) 

end   % plant_dynamics

%% Local helper
function ls = compute_moment_arm(m, params)
% COMPUTE_MOMENT_ARM  Instantaneous CoM-to-nozzle moment arm  (m).
%
% As propellant burns, the CoM shifts upward (for a vehicle with tanks
% located above the nozzle).  The moment arm ls increases accordingly.
%
% Linear model between fully-loaded and dry conditions:
%   frac = (m − ms) / mp          ∈ [0, 1]  (propellant fraction remaining)
%   ls(frac) = ls_dry + (ls_full − ls_dry) · frac
%            → ls_full  when fully loaded  (frac = 1, CoM lower)
%            → ls_dry   when propellant exhausted  (frac = 0, CoM higher)
%
% Inputs:
%   m       (scalar)  current total mass (kg)
%   params            params struct; must contain ms, mp, ls_nominal (m)
%
% Outputs:
%   ls      (scalar)  moment arm from current CoM to nozzle exit (m)


% Propellant fraction remaining — clamp to [0, 1] for robustness
frac = max(0, min(1, (m - params.ms) / params.mp));

% CoM shifts ~5 cm toward the nozzle as propellant burns (geometry estimate).
% ls_full: shorter arm when tanks are full  (CoM is lower, closer to nozzle)
% ls_dry : longer  arm when tanks are empty (CoM is higher, further from nozzle)
ls_full = params.ls_nominal;          % fully loaded    (≈ 0.60 m nominal)
ls_dry  = params.ls_nominal + 0.05;   % propellant gone (CoM 5 cm further)

% Interpolate:  frac=1 → ls_full;  frac=0 → ls_dry
ls = ls_dry + (ls_full - ls_dry) * frac;

end   % compute_moment_arm