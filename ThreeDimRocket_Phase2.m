function [problem,guess] = ThreeDimRocket_Phase2(problem_mp, guess_mp)
% VTVL Rocket (Four-phase Formulation)
% Syntax:  [problem,guess] = myProblem
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%
% Other m-files required: none
% MAT-files required: none
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk
%
% Altered by L. White, University of Sheffield May 2026

% Variables
rho = 1.293;
diameter = 0.15;
CDxy = 1.17;
CDz = 0.61;
gravsl = 9.81;
ls = 1.2;
c = 1488;
rhof = 786;
rhoo = 785;
ms = 8;
mp = 8.88;
hf = 0.77;
ho = 0.28;
ratio = 3.34;
Rearth =6378100;

x_0=0;
x_min=-inf;
x_max = inf;

y_0=0;
y_min=-inf;
y_max = inf;

z_0=0;
z_min=299.5;
z_max=300.5;

vx_0=0;
vx_min=-inf;
vx_max=inf;

vy_0=0;
vy_min=-inf;
vy_max=inf;

vz_0=0;
vz_min=-inf;
vz_max=inf;

m_0=mp+ms;
m_min=ms;
m_max=mp+ms;

% --- TIGHTENED ATTITUDE BOUNDS: Phase 2 hover ---
% Reduced from ±45° to ±8° throughout the hover phase.
% A hovering rocket should remain near-vertical at all times.
% ±8° gives a small buffer above the realistic ~5° hover tolerance
% without being so tight that IPOPT cannot find a feasible solution.
theta_min = -8;
theta_max =  8;
theta_0 = 0;

phi_min   = -8;
phi_max   =  8;
phi_0=0;


p_0=0;
% Angular rate bounds also tightened, a hovering rocket should not
% be spinning. ±5 deg/s gives room for small stabilisation corrections.
p_min = -5;
p_max =  5;

q_0=0;
q_min = -5;
q_max = 5;

T_max=3200;
T_min=T_max*0.01;

alpha_0=0;
alpha_min=-5;
alpha_max=5;

beta_0=0;
beta_min=-5;
beta_max=5;
%------------- BEGIN CODE --------------
% Plant model name, used for Adigator
InternalDynamics=@ThreeDimRocket_Dynamics_Phase2;
SimDynamics=[];

% Analytic derivative files (optional)
problem.analyticDeriv.gradCost=[];
problem.analyticDeriv.hessianLagrangian=[];
problem.analyticDeriv.jacConst=[];

% Settings file
problem.settings=@settings_ThreeDimRocketFourPhase_Phase2;

% Initial time. t0<tf
problem.time.t0_idx=2;
problem.time.t0_min=problem_mp.time.t_min(problem.time.t0_idx);
problem.time.t0_max=problem_mp.time.t_max(problem.time.t0_idx);
guess.t0=guess_mp.time(problem.time.t0_idx);

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_idx=3;
problem.time.tf_min=problem_mp.time.t_min(problem.time.tf_idx);     
problem.time.tf_max=problem_mp.time.t_max(problem.time.tf_idx); 
guess.tf=guess_mp.time(problem.time.tf_idx);

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=problem_mp.parameters.pl;
problem.parameters.pu=problem_mp.parameters.pu;
guess.parameters=guess_mp.parameters;

% Initial conditions for system
problem.states.x0=[];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[x_min y_min z_min vx_min vy_min vz_min m_min theta_min phi_min p_min q_min];
problem.states.x0u=[x_max y_max z_max vx_max vy_max vz_max m_max theta_max phi_max p_max q_max];

% State bounds. xl=< x <=xu
problem.states.xl=[x_min y_min z_min vx_min vy_min vz_min m_min theta_min phi_min p_min q_min];
problem.states.xu=[x_max y_max z_max vx_max vy_max vz_max m_max theta_max phi_max p_max q_max];

% State error bounds
problem.states.xErrorTol_local=[1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2];
problem.states.xErrorTol_integral=[1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2];

% State constraint error bounds
problem.states.xConstraintTol=[1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2 1e-2];


% Terminal state bounds. xfl=< xf <=xfu
% --- FUEL MINIMISATION: Hard terminal position constraint at end of hover ---
pos_tol  = 1.0;
vel_hover = 1.0;   % max exit velocity per axis [m/s] — forces clean MPC entry state

problem.states.xfl = [100-pos_tol, -100-pos_tol, z_min, ...
                       -vel_hover, -vel_hover, -vel_hover, m_min, ...
                       theta_min, phi_min, p_min, q_min];

problem.states.xfu = [100+pos_tol, -100+pos_tol, z_max, ...
                        vel_hover,  vel_hover,  vel_hover, m_max, ...
                       theta_max,   phi_max,    p_max,    q_max];

% Guess the state trajectories with [x0 xf]
guess.states(:,1)=[100 100];
guess.states(:,2)=[-100 -100];
guess.states(:,3)=[300 300];
guess.states(:,4)=[0 0];
guess.states(:,5)=[0 0];
guess.states(:,6)=[0 0];
guess.states(:,7)=[(m_max+m_min)/1.5 (m_max+m_min)/2];
guess.states(:,8)=[0 0];
guess.states(:,9)=[0 0];
guess.states(:,10)=[0 0];
guess.states(:,11)=[0 0];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.

problem.inputs.N=0;

% Input bounds
problem.inputs.ul=[T_min alpha_min beta_min];
problem.inputs.uu=[T_max alpha_max beta_max];

% Bounds on first control action
problem.inputs.u0l=[T_min alpha_min beta_min];
problem.inputs.u0u=[T_max alpha_max beta_max];

% Input constraint error bounds
problem.inputs.uConstraintTol=[0.1 0.1 0.1];

% Guess the input sequences with [u0 uf]
guess.inputs(:,1)=[T_max/2 T_max/2];
guess.inputs(:,2)=[alpha_0 alpha_0];
guess.inputs(:,3)=[beta_0 beta_0];

% Choose the set-points if required
problem.setpoints.states=[];
problem.setpoints.inputs=[];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.ng_eq=0;
problem.constraints.gTol_eq=[];

problem.constraints.gl=[];
problem.constraints.gu=[];
problem.constraints.gTol_neq=[];



% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=[5];
problem.constraints.bu=[5];
problem.constraints.bTol=[0.01];


% store the necessary problem parameters used in the functions
problem.data=problem_mp.data;
problem.data.rho = rho;
problem.data.diameter = diameter;
problem.data.CDxy = CDxy;
problem.data.CDz = CDz;
problem.data.gravsl = gravsl;
problem.data.ls = ls;
problem.data.c = c;
problem.data.rhof = rhof;
problem.data.rhoo = rhoo;
problem.data.ms = ms;
problem.data.hf = hf;
problem.data.ho = ho;
problem.data.ratio = ratio;
problem.data.Rearth = Rearth;

% Get function handles and return to Main.m
problem.data.InternalDynamics=InternalDynamics;
problem.data.functionfg=@fg;
problem.data.plantmodel = func2str(InternalDynamics);
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.sim.functions=SimDynamics;
problem.sim.inputX=[];
problem.sim.inputU=1:length(problem.inputs.ul);
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc,@b_unscaled};
problem.data.functions_unscaled=problem.functions_unscaled;
problem.data.ng_eq=problem.constraints.ng_eq;
problem.constraintErrorTol=[problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];

%------------- END OF CODE --------------

function stageCost=L_unscaled(x,xr,u,ur,p,t,vdat)

% L_unscaled - Returns the stage cost.
% The function must be vectorized and
% xi, ui are column vectors taken as x(:,i) and u(:,i) (i denotes the i-th
% variable)
% 
% Syntax:  stageCost = L(x,xr,u,ur,p,t,data)
%
% Inputs:
%    x  - state vector
%    xr - state reference
%    u  - input
%    ur - input reference
%    p  - parameter
%    t  - time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    stageCost - Scalar or vectorized stage cost
%
%  Remark: If the stagecost does not depend on variables it is necessary to multiply
%          the assigned value by t in order to have right vector dimesion when called for the optimization. 
%          Example: stageCost = 0*t;

%------------- BEGIN CODE --------------

% FUEL MINIMISATION — Lagrange (running) cost
%
% Integrates the propellant mass flow rate dm/dt = T/c [kg/s] over time.
% The integral of this over the full trajectory equals total fuel burned [kg].
%
% This is the standard fuel-optimal objective used by methods such as
% lossless convexification, enabling direct comparison.
%
% Previously: stageCost = 0*t  (no running cost)
% Now:        stageCost = T./c (fuel burn rate)

T = u(:,1);    % thrust magnitude [N]
alpha = u(:,2);
beta  = u(:,3);
c = vdat.c;    % effective exhaust velocity [m/s]  (= 1488 m/s for this rocket)

epsilon_fuel   = 1e-6;
epsilon_gimbal = 1e-2;    % was 1e-4, smooths gimbal chatter significantly

stageCost = T./c + epsilon_fuel.*T.^2 + epsilon_gimbal.*(alpha.^2 + beta.^2);
%------------- END OF CODE --------------


function boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat) 

% E_unscaled - Returns the boundary value cost
%
% Syntax:  boundaryCost=E(x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    boundaryCost - Scalar boundary cost
%
%------------- BEGIN CODE --------------

% BOUNDARY COST — set to zero for pure fuel minimisation
%
% Previously: penalised squared deviation from target position/velocity.
% Now: waypoint accuracy is enforced by hard terminal state bounds
% (set in problem.states.xfl / xfu below), so no soft penalty is needed.
% This gives the same structure as lossless convexification for comparison.

boundaryCost = 0;
%------------- END OF CODE --------------

function bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)

% b_unscaled - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
%
% Syntax:  bc=b(x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
%          
% Output:
%    bc - column vector containing the evaluation of the boundary function 
%
%------------- BEGIN CODE --------------
varargin=varargin{1};
bc=[tf-t0];
%------------- END OF CODE --------------
% When adpative time interval add constraint on time
%------------- BEGIN CODE --------------
if length(varargin)==2
    options=varargin{1};
    t_segment=varargin{2};
    if ((strcmp(options.discretization,'hpLGR')) || (strcmp(options.discretization,'globalLGR')))  && options.adaptseg==1 
        if size(t_segment,1)>size(t_segment,2)
            bc=[bc;diff(t_segment)];
        else
            bc=[bc;diff(t_segment)'];
        end
    end
end

%------------- END OF CODE --------------

