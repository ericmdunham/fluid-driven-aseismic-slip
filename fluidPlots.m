% script to reproduce plots in
% Dunham, E. M., Fluid-driven aseismic fault slip with permeability enhancement and dilatancy, submitted.

savePlots = true; % set to true to save .eps figures

epsilon = [0 1 10 100]; M = length(epsilon);
lambda = logspace(-2,1,10000)'; N = length(lambda);
T2 = zeros(N,M); T2_BV = zeros(N,1);
T3 = zeros(N,M); T2_SL = zeros(N,1);
T3lambdaSmall = zeros(N,M);

for m=1:M
    for n=1:N
        T2(n,m) = fluidDrivenAseismicSlip.evalT2(lambda(n),epsilon(m));
        T3(n,m) = fluidDrivenAseismicSlip.evalT3(lambda(n),epsilon(m));
        if m==1, T2_BV(n,1) = fluidDrivenAseismicSlip.evalT2_BV(lambda(n)); end
        if m==1, T3_SL(n,1) = fluidDrivenAseismicSlip.evalT3_SL(lambda(n)); end
        %if T(n,m)<0, T(n,m)=0; end
        if isnan(T2(n,m)), T2(n,m)=-n; end % need value for interpolation
        if isnan(T3(n,m)), T3(n,m)=-n; end % need value for interpolation
        if isinf(T2(n,m)), T2(n,m)=-n; end % need value for interpolation
        if isinf(T3(n,m)), T3(n,m)=-n; end % need value for interpolation
        T3lambdaSmall(n,m) = 1/lambda(n)^2+1-epsilon(m)-log(4); % asymptotics
    end
end

figure(1),clf
semilogx(lambda,T2,lambda,T2_BV,'--')
xlabel('\lambda = a(t)/(4\alphat)^{1/2}')
ylabel('stress-injection parameter, T')
ylim([0 1])
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'BV19'; legend(str)
title('2D constant pressure injection')
if savePlots,print -depsc2 Tlambda2D, end

figure(2),clf
loglog(lambda,T3,lambda,T3_SL,'--',lambda,T3lambdaSmall,'k--')
xlabel('\lambda = a(t)/(4\alphat)^{1/2}')
ylabel('stress-injection parameter, T')
ylim([1e-1 1e2])
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'SLBV22'; legend(str)
title('3D constant rate injection')
if savePlots,print -depsc2 Tlambda3D, end

% nondimensional pressure
xi = [linspace(0,0.01,10) linspace(0.01,1,100) linspace(1,2,100)]'; I = 1:110;
Tplot2D = 0.5;
for m=1:M
    la2(m) = interp1(T2(:,m),lambda,Tplot2D);
    g2(:,m) = fluidDrivenAseismicSlip.evalG2(xi,la2(m),epsilon(m));
end
g2(I(end)+1,:) = 0;
la2_BV = interp1(T2_BV,lambda,Tplot2D);
g2_BV = fluidDrivenAseismicSlip.evalG2_BV(xi,la2_BV);
g2_BV(I(end)+1) = 0;

Tplot3D = 0.5;
for m=1:M
    la3(m) = interp1(T3(:,m),lambda,Tplot3D);
    g3(:,m) = fluidDrivenAseismicSlip.evalG3(xi,la3(m),epsilon(m));
end
g3(I(end)+1,:) = 0;
la3_SL = interp1(T3_SL,lambda,Tplot3D);
g3_SL = fluidDrivenAseismicSlip.evalG3_SL(xi,la3_SL);
g3_SL(I(end)+1) = 0;

figure(3)
plot(xi,g2,xi,g2_BV,'--')
xlabel('\xi = x / \lambda(4\alphat)^{1/2}')
ylabel('pressure change, p/\Deltap = g(\xi)')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'BV19'; legend(str)
title(['2D constant pressure injection, T = ' num2str(Tplot2D)])
if savePlots,print -depsc2 p2D, end

figure(4)
plot(xi,g3,xi,g3_SL,'--')
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('pressure change, p/\Deltap = g(\xi)')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'SLBV22'; legend(str)
title(['3D constant rate injection, T = ' num2str(Tplot3D)])
ylim([-0.5 2])
if savePlots,print -depsc2 p3D, end

% nondimensional stress drop

figure(5)
plot(xi(I),g2(I,:)-Tplot2D,xi(I),g2_BV(I)-Tplot2D,'--')
line([0 1],[0 0],'color','k','linestyle','--')
xlabel('\xi = x / \lambda(4\alphat)^{1/2}')
ylabel('stress drop, \Delta\tau/(f\Deltap) = g(\xi)-T')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'BV19'; legend(str)
title(['2D constant pressure injection, T = ' num2str(Tplot2D)])
if savePlots,print -depsc2 tau2D, end

figure(6)
plot(xi(I),g3(I,:)-Tplot3D,xi(I),g3_SL(I)-Tplot3D,'--')
line([0 1],[0 0],'color','k','linestyle','--')
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('stress drop, \Delta\tau/(f\Deltap) = g(\xi)-T')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'SLBV22'; legend(str)
title(['3D constant rate injection, T = ' num2str(Tplot3D)])
if savePlots,print -depsc2 tau3D, end

% repeat 3D for larger T

Tplot3D = 16;
for m=1:M
    la3(m) = interp1(T3(:,m),lambda,Tplot3D);
    g3(:,m) = fluidDrivenAseismicSlip.evalG3(xi,la3(m),epsilon(m));
end
g3(I(end)+1,:) = 0;

figure(7)
plot(xi,g3)
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('pressure change, p/\Deltap = g(\xi)')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
title(['3D constant rate injection, T = ' num2str(Tplot3D)])
xlim([0 1.2])
if savePlots,print -depsc2 p3D_largeT, end

% nondimensional stress drop

figure(8)
plot(xi(I),g3(I,:)-Tplot3D)
line([0 1],[0 0],'color','k','linestyle','--')
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('stress drop, \Delta\tau/(f\Deltap) = g(\xi)-T')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
title(['3D constant rate injection, T = ' num2str(Tplot3D)])
if savePlots,print -depsc2 tau3D_largeT, end

% repeat 3D for different T, all epsilon=0

T = [1 3 10 30]; nT = length(T);
for n=1:nT
    la3T(n) = interp1(T3(:,1),lambda,T(n));
    g3T(:,n) = fluidDrivenAseismicSlip.evalG3(xi,la3T(n),0);
end
g3T(I(end)+1,:) = 0;

figure(9)
for n=1:nT
    plot(xi,g3T(:,n)/T(n))
    if n==1,hold on,end
end
hold off
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('pressure change, p/(T\Deltap) = g(\xi)/T')
legend(num2str(T'))
lg=legend(strcat('T = ',num2str(T')));
title('3D constant rate injection, \epsilon = 0')
ylim([0 5]),xlim([0 1.2])
if savePlots,print -depsc2 p3D_varyT, end

%la3_SL = interp1(T3_SL,lambda,Tplot3D);
%g3_SL = fluidDrivenAseismicSlip.evalG3_SL(xi,la3_SL);
%g3_SL(I(end)+1) = 0;

figure(10)
for n=1:nT
    plot([xi*la3T(n); 1.2],[g3T(:,n); 0])
    if n==1,hold on,end
end
plot(xi*la3_SL,g3_SL,'--')
hold off
xlabel('r/(4\alphat)^{1/2}')
ylabel('pressure change, p/\Deltap = g(\xi)')
legend(num2str(T'))
lg=legend(strcat('T = ',num2str(T')));
str = lg.String; str{nT+1} = 'SLBV22'; legend(str)
title('3D constant rate injection, \epsilon = 0')
ylim([0 40]),xlim([0 1.2])
if savePlots,print -depsc2 p3D_varyT3_alt, end

% Cooper Basin estimates for 2012 stimulation of Habanero 4 well

CB.Q = 0.025; % 25 L/s
CB.k = 2e-13;
CB.eta = 4e-4;
CB.w = 6;
CB.Deltap = CB.Q*CB.eta/(4*pi*CB.k*CB.w)*1e-6;
CB.tau = 10.3;
CB.sigmaprime = 28;
CB.f = 0.6;
CB.T = (CB.f*CB.sigmaprime-CB.tau)/(CB.f*CB.Deltap);
disp('Cooper Basin:')
disp(['Deltap = ' num2str(CB.Deltap) ' MPa'])
disp(['T = ' num2str(CB.T)])
