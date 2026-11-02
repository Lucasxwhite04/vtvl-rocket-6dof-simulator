function params = get_params()
% GET_PARAMS  Returns params struct inside Simulink MATLAB Function blocks.
% Simulink MATLAB Function blocks cannot read the workspace directly.
% This function acts as a bridge, params is loaded here from the base
% workspace using evalin.
%
% run_nominal.m must have called plant_params() and assigned
% the result to 'params' in the base workspace before sim() is called.
% L. White, University of Sheffield, May 2026
params = evalin('base', 'params');
end