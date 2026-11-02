function [problem,guess,phaseoptions] = ThreeDimRocketFourPhase
% VTVL Rocket (Four-phase Formulation)
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

% Initial and final time for different phases. Let t_min(end)=t_max(end) if tf is fixed.
% 4 Phase - With Hover
problem.mp.time.t_min=[0 0 0 0 0];     
problem.mp.time.t_max=[0 20 30 50 60 ]; 
guess.mp.time=[0 10 15 35 40];
% % New: 4 time nodes for 3 phases (Ascent → Descent → Land)
% problem.mp.time.t_min = [0   0   0   0];
% problem.mp.time.t_max = [0  20  40  60];
% guess.mp.time         = [0  10  35  42];


% Parameters bounds. pl=< p <=pu
problem.mp.parameters.pl=[];
problem.mp.parameters.pu=[];
guess.mp.parameters=[];

% Bounds for linkage boundary constraints bll =< bclink(x0,xf,u0,uf,p,t0,tf,vdat) =< blu
% 4 Phase Hover
% Linkage: 4 phases = 3 linkages × (11 states + 3 inputs + 1 time) = 45 entries
problem.mp.constraints.bll.linear = zeros(1, 45);
problem.mp.constraints.blu.linear = zeros(1, 45);
problem.mp.constraints.blTol.linear = 0.01 * ones(1, 45);
% % 3 Phase Without Hover
% % Linkage now has 11 states + 3 inputs + 1 time = 15 components per linkage,
% % times 2 linkages = 30 entries
% problem.mp.constraints.bll.linear = zeros(1, 30);
% problem.mp.constraints.blu.linear = zeros(1, 30);
% problem.mp.constraints.blTol.linear = 0.01 * ones(1, 30);

problem.mp.constraints.bll.nonlinear=[];
problem.mp.constraints.blu.nonlinear=[];
problem.mp.constraints.blTol.nonlinear=[];

% Get function handles
problem.mp.linkfunctions=@bclink;

% store the necessary problem parameters used in the functions
problem.mp.data = [];

% define different phases of OCP
% 4 Phase - With Hover
[problem.phases{1},guess.phases{1}] = ThreeDimRocket_Phase1(problem.mp, guess.mp);
[problem.phases{2},guess.phases{2}] = ThreeDimRocket_Phase2(problem.mp, guess.mp);
[problem.phases{3},guess.phases{3}] = ThreeDimRocket_Phase3(problem.mp, guess.mp);
[problem.phases{4},guess.phases{4}] = ThreeDimRocket_Phase4(problem.mp, guess.mp);
phaseoptions{1}=settings_ThreeDimRocketFourPhase_Phase1(20);
phaseoptions{2}=settings_ThreeDimRocketFourPhase_Phase2(15);
phaseoptions{3}=settings_ThreeDimRocketFourPhase_Phase3(20);
phaseoptions{4}=settings_ThreeDimRocketFourPhase_Phase4(20);
% % Drop the Phase 2 (hover) instantiation
% [problem.phases{1},guess.phases{1}] = ThreeDimRocket_Phase1(problem.mp, guess.mp);
% [problem.phases{2},guess.phases{2}] = ThreeDimRocket_Phase3(problem.mp, guess.mp);  % was Phase 3
% [problem.phases{3},guess.phases{3}] = ThreeDimRocket_Phase4(problem.mp, guess.mp);  % was Phase 4
% phaseoptions{1} = settings_ThreeDimRocketFourPhase_Phase1(40);   % more nodes for ascent
% phaseoptions{2} = settings_ThreeDimRocketFourPhase_Phase3(50);   % most nodes for descent
% phaseoptions{3} = settings_ThreeDimRocketFourPhase_Phase4(40);   % nodes for landing burn
%------------- END OF CODE --------------


function [blc_linear, blc_nonlinear]=bclink(x0,xf,u0,uf,p,t0,tf,vdat)

% bclink - Returns the evaluation of the linkage boundary constraints: bll =< bclink(x0,xf,u0,uf,p,t0,tf,vdat) =< blu
%
% Syntax:  [blc_linear, blc_nonlinear]=bclink(x0,xf,u0,uf,p,t0,tf,vdat)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    vdat- structured variable containing the values of additional data used inside
%          the function
%
%          
% Output:
%    blc_linear - column vector containing the evaluation of the linear linkage boundary constraint functions
%    blc_nonlinear - column vector containing the evaluation of the nonlinear linkage boundary constraint functions
%
%------------- BEGIN CODE --------------

% 4 Phase
blc_linear=[xf{1}(1)-x0{2}(1);xf{1}(2)-x0{2}(2);xf{1}(3)-x0{2}(3);xf{1}(4)-x0{2}(4);xf{1}(5)-x0{2}(5);xf{1}(6)-x0{2}(6);xf{1}(7)-x0{2}(7);xf{1}(8)-x0{2}(8);xf{1}(9)-x0{2}(9);xf{1}(10)-x0{2}(10);xf{1}(11)-x0{2}(11);xf{2}(1)-x0{3}(1);xf{2}(2)-x0{3}(2);xf{2}(3)-x0{3}(3);xf{2}(4)-x0{3}(4);xf{2}(5)-x0{3}(5);xf{2}(6)-x0{3}(6);xf{2}(7)-x0{3}(7);xf{2}(8)-x0{3}(8);xf{2}(9)-x0{3}(9);xf{2}(10)-x0{3}(10);xf{2}(11)-x0{3}(11);xf{3}(1)-x0{4}(1);xf{3}(2)-x0{4}(2);xf{3}(3)-x0{4}(3);xf{3}(4)-x0{4}(4);xf{3}(5)-x0{4}(5);xf{3}(6)-x0{4}(6);xf{3}(7)-x0{4}(7);xf{3}(8)-x0{4}(8);xf{3}(9)-x0{4}(9);xf{3}(10)-x0{4}(10);xf{3}(11)-x0{4}(11);uf{1}(1)-u0{2}(1);uf{1}(2)-u0{2}(2);uf{1}(3)-u0{2}(3);uf{2}(1)-u0{3}(1);uf{2}(2)-u0{3}(2);uf{2}(3)-u0{3}(3);uf{3}(1)-u0{4}(1);uf{3}(2)-u0{4}(2);uf{3}(3)-u0{4}(3);tf(1)-t0(2);tf(2)-t0(3);tf(3)-t0(4)];
% New 3 Phase - No hover
% In bclink: drop the third phase linkage block
% blc_linear = [...
%     xf{1}(1)-x0{2}(1); xf{1}(2)-x0{2}(2); xf{1}(3)-x0{2}(3); ...
%     xf{1}(4)-x0{2}(4); xf{1}(5)-x0{2}(5); xf{1}(6)-x0{2}(6); ...
%     xf{1}(7)-x0{2}(7); xf{1}(8)-x0{2}(8); xf{1}(9)-x0{2}(9); ...
%     xf{1}(10)-x0{2}(10); xf{1}(11)-x0{2}(11); ...
%     xf{2}(1)-x0{3}(1); xf{2}(2)-x0{3}(2); xf{2}(3)-x0{3}(3); ...
%     xf{2}(4)-x0{3}(4); xf{2}(5)-x0{3}(5); xf{2}(6)-x0{3}(6); ...
%     xf{2}(7)-x0{3}(7); xf{2}(8)-x0{3}(8); xf{2}(9)-x0{3}(9); ...
%     xf{2}(10)-x0{3}(10); xf{2}(11)-x0{3}(11); ...
%     uf{1}(1)-u0{2}(1); uf{1}(2)-u0{2}(2); uf{1}(3)-u0{2}(3); ...
%     uf{2}(1)-u0{3}(1); uf{2}(2)-u0{3}(2); uf{2}(3)-u0{3}(3); ...
%     tf(1)-t0(2); tf(2)-t0(3)];

blc_nonlinear=[];
%------------- END OF CODE --------------
