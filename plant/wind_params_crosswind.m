function p = wind_params_crosswind(speed_mps, direction_deg)
% WIND_PARAMS_CROSSWIND  Wind parameters for the crosswind disturbance scenario.
%
% Constant mean wind plus moderate turbulence.
%
% Inputs:
%   speed_mps     mean wind speed at 10 m reference altitude  (m/s)
%       
%   direction_deg direction wind blows FROM, meteorological convention  (deg)
%                 e.g. 270 = from the west, blowing in the +X direction
%
% Outputs:
%   p  struct  wind parameter struct accepted by wind_model.m
%
% Usage:
%   p = wind_params_crosswind(5, 270);
%   w = wind_model(t, z, p, 42);
% L. White, University of Sheffield, May 2026


% Assigns the user-provided speed to the mean_speed parameter.
p.mean_speed = speed_mps;    

% Assigns the user-provided wind direction.
p.direction_deg = direction_deg; 

% Sets the reference altitude for the wind shear power law to 10 meters.
p.ref_altitude          = 10;           

% Sets the Hellmann exponent for the wind shear profile.
p.shear_exp             = 0.14;         

% Sets the turbulence scaling multiplier to 0.15 (15%).
p.turbulence_intensity  = 0.15;         

% Sets the standard deviation (base strength) of the gusts to 0.5 m/s.
p.gust_sigma            = 0.5;          

% Sets the frequency spectrum width for the gusts to 0.5 rad/s.
p.gust_bandwidth        = 0.5;          

end