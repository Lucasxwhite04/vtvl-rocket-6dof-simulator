function z_meas = baro_model(x_true, params, t, seed)
% BARO_MODEL  Simulated barometric altimeter.
%
% Altitude-only measurement at ~50 Hz with a fixed turn-on bias and
% per-sample white noise.
% L. White, University of Sheffield, May 2026

%%  PERSISTENT BIAS - generated once per seed, not every timestep
persistent cached_seed bias_rt

p_baro = params.baro;

% Check if a new simulation run (with a new seed) has started.
if isempty(cached_seed) || cached_seed ~= seed
    
    % Seed the random number generator uniquely for the barometer.
    rng(seed + 7000, 'twister');   
    
    % Calculate the permanent altitude offset for this specific flight.
    bias_rt     = p_baro.bias + p_baro.bias_sigma * randn();
    
    cached_seed = seed;
end

%%  MEASUREMENT
% Extract the true altitude from the 14-state truth vector.
z_true = x_true(3);

% Generate the high-frequency measurement noise.
noise = p_baro.noise_std * randn();

% Construct the final corrupted sensor reading.
z_meas = z_true + bias_rt + noise;   % (scalar) m

end   % baro_model