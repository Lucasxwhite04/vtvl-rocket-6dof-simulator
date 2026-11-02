function w = wind_model(t, z, params_wind, seed)
% WIND_MODEL  Wind velocity profile in the inertial frame.
%
% Produces a deterministic (seeded), repeatable wind signal composed of:
%   1. Mean wind shear:  power-law altitude profile
%   2. Turbulence:       Dryden-like band-limited noise via sum-of-sinusoids
%                        (Shinozuka method — deterministic given seed)
% L. White, University of Sheffield, May 2026


%  PERSISTENT TURBULENCE STATE
%'persistent' variables retain their values in memory between function 
% calls, acting like private global variables for this function.
persistent cached_seed cached_bandwidth cached_sigma ...
           phases_x phases_y phases_z omegas amplitudes

% Fills in missing fields with safe default values.
params_wind = apply_defaults(params_wind);

% --- Decide whether to (re)generate turbulence phases ---
% Checks if this is the very first time the function is called (isempty), 
% or if the Monte Carlo runner has changed the seed or turbulence parameters.
need_reinit = isempty(cached_seed)                        || ...
              cached_seed      ~= seed                    || ...
              cached_bandwidth ~= params_wind.gust_bandwidth || ...
              cached_sigma     ~= params_wind.gust_sigma;

if need_reinit
    % Set up a grid of 20 frequencies (omegas) and constant amplitudes.
    N_harmonics = 20;
    omega_base = params_wind.gust_bandwidth;    % rad/s
    omega_max  = 10 * omega_base;               % rad/s
    omegas     = linspace(0.1 * omega_base, omega_max, N_harmonics);  
    
    % Calculate the amplitude for each sine wave based on desired gust intensity.
    amplitudes = params_wind.gust_sigma * sqrt(2 / N_harmonics) * ones(1, N_harmonics);
    
    % Seed the random number generator and generate random phase shifts.
    rng(seed, 'twister');
    phases_x = 2*pi * rand(1, N_harmonics);   
    phases_y = 2*pi * rand(1, N_harmonics);
    phases_z = 2*pi * rand(1, N_harmonics);
    
    % Save the current setup to the persistent variables so we don't 
    %       trigger 'need_reinit' on the next timestep.
    cached_seed      = seed;
    cached_bandwidth = params_wind.gust_bandwidth;
    cached_sigma     = params_wind.gust_sigma;
end

%%  MEAN WIND SHEAR  (power-law altitude profile)
% Clamp the minimum altitude to 1.0 meters.
z_eff = max(z, 1.0);

% Calculate the average wind speed at the current altitude using a power law.
if params_wind.mean_speed > 0
    U_mean = params_wind.mean_speed * ...
             (z_eff / params_wind.ref_altitude) ^ params_wind.shear_exp;
else
    U_mean = 0;
end

% Convert meteorological wind direction (degrees) into X/Y vectors.
dir_rad = params_wind.direction_deg * pi / 180;
wx_mean = -U_mean * cos(dir_rad);   
wy_mean = -U_mean * sin(dir_rad);   
wz_mean = 0;                        

%%  TURBULENCE  (Dryden-like, sum-of-sinusoids)
% Evaluate the 20 sine waves for the *current time* (t) and sum them up.
wx_turb = sum(amplitudes .* cos(omegas * t + phases_x));
wy_turb = sum(amplitudes .* cos(omegas * t + phases_y));

% Vertical gusts are halved (* 0.5).
wz_turb = sum(amplitudes .* cos(omegas * t + phases_z)) * 0.5;

% Create an altitude multiplier that goes from 0.0 at ground to 1.0 at 50m.
alt_factor = min(z_eff / 50.0, 1.0);
intensity  = params_wind.turbulence_intensity * alt_factor;

% Apply the intensity and altitude damping to the raw turbulent signal
wx_turb = intensity * wx_turb;
wy_turb = intensity * wy_turb;
wz_turb = intensity * wz_turb;

%%  ASSEMBLE OUTPUT
% Add the static mean wind (shear) to the varying gusts (turbulence).
w = [ wx_mean + wx_turb ;
      wy_mean + wy_turb ;
      wz_mean + wz_turb ];
end   

%%  LOCAL HELPER
function p = apply_defaults(p) %% THIS IS A BACKUP
    % Checks if a field exists, and if not, creates it with a default value.
    if ~isfield(p, 'mean_speed'),           p.mean_speed           = 0;     end  
    if ~isfield(p, 'direction_deg'),        p.direction_deg        = 270;   end  
    if ~isfield(p, 'ref_altitude'),         p.ref_altitude         = 10;    end  
    if ~isfield(p, 'shear_exp'),            p.shear_exp            = 0.14;  end  
    if ~isfield(p, 'turbulence_intensity'), p.turbulence_intensity = 0.0;   end  
    if ~isfield(p, 'gust_sigma'),           p.gust_sigma           = 0.0;   end  
    if ~isfield(p, 'gust_bandwidth'),       p.gust_bandwidth       = 0.5;   end  
end