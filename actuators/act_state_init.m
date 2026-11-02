function act_in = act_state_init()
% ACT_STATE_INIT  Initialise TVC servo persistent state.
% Call once at simulation startup before first tvc_servo call.
% All states start at zero — gimbal centred, at rest.
% L. White, University of Sheffield, May 2026

act_in.s_alpha      = [0; 0];   % pitch axis state [angle; rate] (rad, rad/s)
act_in.s_beta       = [0; 0];   % roll  axis state [angle; rate] (rad, rad/s)
act_in.alpha_cmd_rl = 0;        % last rate-limited pitch command (rad)
act_in.beta_cmd_rl  = 0;        % last rate-limited roll  command (rad)
act_in.alpha_prev   = 0;        % last pitch output (rad)
act_in.beta_prev    = 0;        % last roll  output (rad)
end