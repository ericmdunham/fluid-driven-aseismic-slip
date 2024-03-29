% script to reproduce estimates of T, lambda, etc., for 2012 stimulation of Habanero 4
% Dunham, E. M., Fluid-driven aseismic fault slip with permeability enhancement and dilatancy, submitted.

eta = 8e-4; % Pa s
k = 1e-12; % m^2
Q = 25e-3; % m^3/s
S = 1e-2*5e-8; % 1/Pa
alpha = k/(S*eta); % m^2/s
w = 6; % m

Deltap = Q*eta/(4*pi*k*w)*1e-6 % MPa

epsilon = 0;

Sv = 97.9; SHmax = 150.9; Shmin = 124.4; p = 73.1; % Barton et al. (2013)
dip = 10; % deg
[stt,snt,snn]=rotate_stress(SHmax,0,Sv,-deg2rad(dip))
sigma = snn; tau = snt; % MPa
sigmaEff = sigma-p; % MPa
f = 0.5;

T = (f*sigmaEff-tau)/(f*Deltap)

lambda = fzero(@(lambda) T-fluidDrivenAseismicSlip.evalT3(lambda,epsilon),0.25)

R = 0.178/2; % m

tD = linspace(1,7); % time (day)
t = tD*60*60*24; % time (s)
xi = R./(lambda*sqrt(4*alpha*t)); % dimensionless argument
gWB = fluidDrivenAseismicSlip.evalG3(xi,lambda,epsilon);
gamma = 0.5772;
gWBapprox = 1/lambda^2+log(lambda^2)-1-epsilon+gamma+4*pi;


figure(1),clf
subplot(2,1,1)
plot(tD,Deltap*gWB,tD,Deltap*gWBapprox+0*tD,'--')
xlabel('time (day)'),ylabel('injector pressure change (MPa)')

subplot(2,1,2)
rSlip = lambda*sqrt(4*alpha*t); % m (700 m after 7 d)
plot(tD,rSlip)
xlabel('time (day)'),ylabel('crack radius (m)')