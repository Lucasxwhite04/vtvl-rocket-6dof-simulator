function p = wind_params_off()
% WIND_PARAMS_OFF Wind parameters for the zero-wind (nominal) scenario.
%
% Returns a params_wind struct with all wind disabled.
% Use for Scenario 1 (nominal — no dispersions, no wind).
%
% Usage:
%   w = wind_model(t, z, wind_params_off(), 0);
%
% Outputs:
%   p struct wind parameter struct accepted by wind_model.m
% L. White, University of Sheffield, May 2026

% Sets the base wind speed to exactly zero.
p.mean_speed = 0;             

% Defines the direction the wind blows from (meteorological convention).
p.direction_deg = 270;        

% Sets the reference altitude for the wind shear power law.
p.ref_altitude = 10;          

% Sets the Hellmann exponent for wind shear profile.
p.shear_exp = 0.14;           

%Sets the turbulence scaling multiplier to zero.
p.turbulence_intensity = 0.0; 

% Sets the standard deviation (strength) of the gusts to zero.
p.gust_sigma = 0.0;           

% Sets the frequency spectrum width for the (disabled) gusts.
p.gust_bandwidth = 0.5;       

end