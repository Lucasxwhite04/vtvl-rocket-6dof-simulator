function [x_hat_new, P_new] = ekf_predict(x_hat, P, accel_meas, gyro_meas, dt, params, T_cmd)
x_hat      = x_hat(:);
accel_meas = accel_meas(:);
gyro_meas  = gyro_meas(:);
% EKF_PREDICT  EKF time update (predict step) driven by IMU measurements.
%
% Propagates the 17-state estimate forward by dt seconds using strapdown
% inertial navigation equations plus rocket mass depletion. Runs at IMU
% rate (200 Hz => dt = 0.005 s).
% L. White, University of Sheffield, May 2026

% --- Unpack state ---
r   = x_hat(1:3);
v   = x_hat(4:6);
q   = x_hat(7:10);    % [qw; qx; qy; qz]
b_g = x_hat(11:13);
b_a = x_hat(14:16);
m   = x_hat(17);      % vehicle mass (kg)

% --- Bias-corrected IMU measurements ---
% Subtract the EKF's current estimate of the sensor biases from the raw data.
omega = gyro_meas  - b_g;   % corrected angular rate (rad/s)
f     = accel_meas - b_a;   % corrected specific force, body frame (m/s^2)

% --- Rotation matrix from quaternion (body -> inertial) ---
% Convert the 4D attitude quaternion into a 3x3 Direction Cosine Matrix.
R = quat2rotm_local(q);   % (3x3)

% --- Strapdown propagation ---
g_vec = [0; 0; -params.g0];

% Rotate the cleaned accelerometer data into the Earth frame.
f_inertial = R * f;

% Add gravity to the specific force to get true coordinate acceleration, 
% then integrate using standard Newtonian kinematics.
a = f_inertial + g_vec;
r_new = r + v*dt + 0.5*a*dt^2;
v_new = v + a*dt;

% Update the rocket's rotation (attitude) over the tiny timestep.
angle = norm(omega) * dt;
if angle > 1e-10
    axis  = omega / norm(omega);
    dq    = [cos(angle/2); sin(angle/2) * axis];
else
    % A mathematical shortcut if the rocket is barely spinning.
    dq = [1; 0.5 * omega * dt];
    dq = dq / norm(dq);
end

% Apply the small rotation to the current rotation via Hamilton product.
q_new = quat_multiply(q, dq);
q_new = q_new / norm(q_new);   % renormalise to prevent numerical drift

% Propagate the biases.
b_g_new = b_g;
b_a_new = b_a;

% --- Mass depletion ---
% Subtract the burned fuel mass using the commanded thrust.
T_mag   = max(T_cmd, 0);            
m_new   = m - (T_mag / params.c) * dt;
m_new   = max(m_new, params.ms);    % cannot fall below dry mass



% --- Assemble new state (17x1) ---
x_hat_new = [r_new; v_new; q_new; b_g_new; b_a_new; m_new];

% --- Linearised state transition matrix F (17x17) ---
% The Jacobian matrix (partial derivatives of the physics equations).
F = eye(17);

% Velocity error integrates into position error
F(1:3, 4:6) = eye(3) * dt;

% Attitude error integrates into velocity error (via rotated specific force)
dRf_dq = dRotVec_dQuat(q, f);   % (3x4)
F(4:6, 7:10) = dRf_dq * dt;

% Accelerometer bias error directly corrupts velocity
F(4:6, 14:16) = -R * dt;

% Gyro bias error directly corrupts attitude (quaternion)
Xi = quat_Xi(q);   % (4x3)
F(7:10, 11:13) = -0.5 * Xi * dt;

% The mass row remains uncoupled (zeros everywhere except the diagonal).

% --- Process noise matrix Q (17x17) ---
% Build the matrix that injects new uncertainty every timestep.
Q = build_process_noise(params, dt);

% --- Covariance propagation ---
% The fundamental equation of the Kalman Filter prediction step.
P_new = F * P * F' + Q;

% A mathematical safeguard to force the matrix to be perfectly symmetrical.
P_new = (P_new + P_new') / 2;
end

%% LOCAL HELPERS
function R = quat2rotm_local(q)
% QUAT2ROTM_LOCAL  Rotation matrix from unit quaternion (body -> inertial).
% q = [qw; qx; qy; qz], scalar-first convention.
qw = q(1); qx = q(2); qy = q(3); qz = q(4);
R = [1-2*(qy^2+qz^2),   2*(qx*qy-qw*qz),   2*(qx*qz+qw*qy);
     2*(qx*qy+qw*qz),   1-2*(qx^2+qz^2),   2*(qy*qz-qw*qx);
     2*(qx*qz-qw*qy),   2*(qy*qz+qw*qx),   1-2*(qx^2+qy^2)];
end

function q_out = quat_multiply(q1, q2)
% QUAT_MULTIPLY  Hamilton product q1 ⊗ q2, scalar-first convention.
w1=q1(1); x1=q1(2); y1=q1(3); z1=q1(4);
w2=q2(1); x2=q2(2); y2=q2(3); z2=q2(4);
q_out = [w1*w2 - x1*x2 - y1*y2 - z1*z2;
         w1*x2 + x1*w2 + y1*z2 - z1*y2;
         w1*y2 - x1*z2 + y1*w2 + z1*x2;
         w1*z2 + x1*y2 - y1*x2 + z1*w2];
end

function Xi = quat_Xi(q)
% QUAT_XI  The Xi matrix such that dq/dt = 0.5 * Xi(q) * omega.
% Maps body angular rate omega (3x1) to quaternion rate (4x1).
% q = [qw; qx; qy; qz]
qw=q(1); qx=q(2); qy=q(3); qz=q(4);
Xi = [-qx, -qy, -qz;
       qw, -qz,  qy;
       qz,  qw, -qx;
      -qy,  qx,  qw];
end

function dRf_dq = dRotVec_dQuat(q, f)
% DROTVEC_DQUAT  Jacobian of R(q)*f with respect to q (3x4).
% Used in the F matrix for the velocity row.
qw=q(1); qx=q(2); qy=q(3); qz=q(4);
fx=f(1); fy=f(2); fz=f(3);
dR_dqw = 2*[ 0,    -qz,   qy;
              qz,   0,   -qx;
             -qy,   qx,   0 ];
dR_dqx = 2*[ 0,    qy,   qz;
              qy, -2*qx, -qw;
              qz,   qw, -2*qx];
dR_dqy = 2*[-2*qy, qx,   qw;
              qx,   0,   qz;
             -qw,   qz, -2*qy];
dR_dqz = 2*[-2*qz, -qw,  qx;
              qw,  -2*qz, qy;
              qx,   qy,   0 ];
dRf_dq = [dR_dqw*[fx;fy;fz], dR_dqx*[fx;fy;fz], ...
          dR_dqy*[fx;fy;fz], dR_dqz*[fx;fy;fz]];
end

function Q = build_process_noise(params, dt)
% BUILD_PROCESS_NOISE  Discrete-time process noise covariance (17x17).
p = params.ekf;
% Continuous-time noise spectral densities
q_accel = p.accel_noise_std.^2;     % (m/s^2)^2 / Hz
q_gyro  = p.gyro_noise_std.^2;      % (rad/s)^2 / Hz
q_bg    = p.gyro_bias_rw.^2;        % gyro bias random walk, (rad/s)^2/s
q_ba    = p.accel_bias_rw.^2;       % accel bias random walk, (m/s^2)^2/s
q_m     = p.mass_process_noise_std^2; % mass process noise, (kg)^2/s

Q = zeros(17);
% Position noise from double-integrating accel noise
Q(1:3, 1:3) = diag(q_accel) * dt^3 / 3;
Q(1:3, 4:6) = diag(q_accel) * dt^2 / 2;
Q(4:6, 1:3) = diag(q_accel) * dt^2 / 2;
% Velocity noise from integrating accel noise
Q(4:6, 4:6) = diag(q_accel) * dt;
% Attitude noise from gyro noise
Q(7:10, 7:10) = diag([q_gyro; q_gyro(1)]) * dt * 0.25;
% Bias random walks
Q(11:13, 11:13) = diag(q_bg) * dt;
Q(14:16, 14:16) = diag(q_ba) * dt;
% Mass process noise (thrust-uncertainty driven)
Q(17, 17) = q_m * dt;
end