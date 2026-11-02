function [deq] = ThreeDimRocket_Dynamics_Phase3(x,u,p,t,vdat)
% VTVL Rocket (Four-phase Formulation)
%
% Syntax:  
%          [dx] = myProblem_Dynamics_Sim(x,u,p,t,vdat)
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    vdat - structured variable containing the values of additional data used inside
%          the function%      
% Output:
%    dx - time derivative of x
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
%------------- BEGIN CODE --------------
x = full(x);
u = full(u);
p = full(p);
% Variables
rho = vdat.rho;
diameter = vdat.diameter;
CDxy = vdat.CDxy;
CDz = vdat.CDz;
gravsl = vdat.gravsl;
ls = vdat.ls;
c = vdat.c;
rhof = vdat.rhof;
rhoo = vdat.rhoo;
ms = vdat.ms;
hf = vdat.hf;
ho = vdat.ho;
ratio = vdat.ratio;
Rearth = vdat.Rearth;

Axy = diameter.*ls;
Az = 3.1415.*(diameter/2).^2;
radius = diameter./2;

% States
X = x(:,1);
Y = x(:,2);
Z = x(:,3);
vex = x(:,4);
vey = x(:,5);
vez = x(:,6);
m = x(:,7);
theta = x(:,8);
phi = x(:,9);
psi = 0;
p = x(:,10);
q = x(:,11);

% Inputs
T = u(:,1);
alpha = u(:,2);
beta = u(:,3);

% Centre of mass
mp = m-ms;

mf = mp./(1+ratio);
mo = mp-mf;

zCoMs = ls/2;

zCoMf = 0.5.*mf./(rhof.*Az)+hf;
zCoMo = 0.5.*mo./(rhoo.*Az)+ho;

zCoM = (zCoMs.*ms+zCoMf.*mf+zCoMo.*mo)./(m);

% Velocity
dx=vex;
dy=vey;
dz=vez;

vrx = vex.*cosd(theta).*cosd(psi)+ vey.*cosd(theta).*sind(psi)+vez.*-sind(theta);
vry = vex.*(sind(phi).*sind(theta).*cosd(psi)-cosd(phi).*sind(psi))+vey.*(sind(phi).*sind(theta).*sind(psi)+cosd(phi).*cosd(psi))+vez.*sind(phi).*cosd(theta);
vrz = vex.*(cosd(phi).*sind(theta).*cosd(psi)+sind(phi).*sind(psi))+vey.*(cosd(phi).*sind(theta).*sind(psi)-sind(phi).*cosd(psi))+vez.*cosd(phi).*cosd(theta);


% New: smooth v·|v| form, C^1 everywhere
DRX = -0.5 .* rho .* vrx .* abs(vrx) .* CDxy .* Axy;
DRY = -0.5 .* rho .* vry .* abs(vry) .* CDxy .* Axy;
DRZ = -0.5 .* rho .* vrz .* abs(vrz) .* CDz  .* Az;

Dex=DRX.*cosd(theta).*cosd(psi)+DRY.*(sind(phi).*sind(theta).*cosd(psi)-cosd(phi).*sind(psi))+DRZ.*(cosd(phi).*sind(theta).*cosd(psi)+sind(phi).*sind(psi));
Dey=DRX.*cosd(theta).*sind(psi)+DRY.*(sind(phi).*sind(theta).*sind(psi)+cosd(phi).*cosd(psi))+DRZ.*(cosd(phi).*sind(theta).*sind(psi)-sind(phi).*cosd(psi));
Dez=DRX.*-sind(theta)+DRY.*sind(phi).*cosd(theta)+DRZ.*cosd(phi).*cosd(theta);

% Thrust
Trx = T.*sind(alpha);
Try = T.*sind(beta);
Trz = T.*cosd(alpha).*cosd(beta);

Tex=-Trx.*cosd(theta).*cosd(psi)+Try.*(sind(phi).*sind(theta).*cosd(psi)-cosd(phi).*sind(psi))+Trz.*(cosd(phi).*sind(theta).*cosd(psi)+sind(phi).*sind(psi));
Tey=-Trx.*cosd(theta).*sind(psi)+Try.*(sind(phi).*sind(theta).*sind(psi)+cosd(phi).*cosd(psi))+Trz.*(cosd(phi).*sind(theta).*sind(psi)-sind(phi).*cosd(psi));
Tez=-Trx.*-sind(theta)+Try.*sind(phi).*cosd(theta)+Trz.*cosd(phi).*cosd(theta);

% Gravity
grav = gravsl.*(1+Z/Rearth).^-2;

%Acceleration
dvx = (1./m).*(Tex+Dex);
dvy = (1./m).*(Tey+Dey);
dvz = (1./m).*(Tez+Dez)-grav;

% Mass flow rate
dm = -T./c;

% Inertia
lf = mf./(rhof.*Az);
lo = mo./(rhoo.*Az);

Ixy = (1/4).*(m.*(radius)^2)+(1/12).*(ms.*ls.^2+mf.*lf.^2+mo.*lo.^2)+ms.*(zCoMs-zCoM).^2+mf.*(zCoMf-zCoM).^2+mo.*(zCoMo-zCoM).^2;

% Angular velocity
dtheta = q;
dphi = p;

% Angular acceleration
dp = Try.*zCoM./Ixy;
dq = Trx.*zCoM./Ixy; 


% Differential
deq=[dx,dy,dz,dvx,dvy,dvz,dm,dtheta,dphi,dp,dq];

%------------- END OF CODE --------------