function [r_meas, v_meas, valid] = gps_model(x_true, params, t, seed)
% GPS_MODEL  Simulated GPS: position + velocity in inertial frame.
%
% Runs at low rate (10 Hz default). Returns valid = false between updates
% so the EKF skips the measurement update and coasts on IMU prediction.
% L. White, University of Sheffield, May 2026

%%  PERSISTENT DROPOUT SCHEDULE - generated once per seed
persistent cached_seed dropout_schedule n_epochs_generated
p_gps = params.gps;

% Determine exactly how many GPS samples will occur during the mission.
mission_dur = p_gps.mission_duration_s;   % s, e.g. 60
n_epochs    = ceil(mission_dur * p_gps.rate_hz) + 10;   % small margin

% Check if we need to regenerate the dropout schedule.
need_reinit = isempty(cached_seed) || cached_seed ~= seed || ...
              n_epochs > n_epochs_generated;

if need_reinit
    % Seed the RNG uniquely (+5000) for the GPS system.
    rng(seed + 5000, 'twister');   
    
    % Create an array of True/False values for every GPS epoch.
    dropout_schedule    = rand(1, n_epochs) < p_gps.dropout_prob;
    
    n_epochs_generated  = n_epochs;
    cached_seed         = seed;
end

%%  CHECK WHETHER THIS TIMESTEP HAS A GPS UPDATE
% Check if the current time 't' aligns with a GPS update window.
gps_period     = 1 / p_gps.rate_hz;                        % s
is_update_step = mod(t, gps_period) < (0.5 / 1000 + 1e-9); % within 0.5 ms of epoch

% Determine which specific GPS epoch we are currently in.
epoch_idx = min(floor(t * p_gps.rate_hz) + 1, n_epochs_generated);
is_dropped = dropout_schedule(epoch_idx);

% The measurement is only valid if we are on an update step AND we didn't drop out.
valid = is_update_step && ~is_dropped;

% Return empty zeros if the GPS is not valid this timestep.

if ~valid
    r_meas = zeros(3, 1);
    v_meas = zeros(3, 1);
    return;
end

%%  NOISY MEASUREMENT
% Add white noise to the true physical position and velocity.
r_meas = x_true(1:3) + p_gps.pos_noise_std .* randn(3, 1);   % (3×1) m
v_meas = x_true(4:6) + p_gps.vel_noise_std .* randn(3, 1);   % (3×1) m/s

end   % gps_model